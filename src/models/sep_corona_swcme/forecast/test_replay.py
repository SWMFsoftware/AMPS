#!/usr/bin/python3
"""Independent image-space and linked-runtime checks for shock replay."""
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

import numpy as np
from scipy.integrate import solve_ivp

import replay
from manufactured_event import AMBIENT, REFERENCE, TRUTH, synthetic, timestamp

CORE = Path(sys.argv.pop(1)).resolve()


class ForecastTest(unittest.TestCase):
    def test_drag_arrival_against_independent_ode(self):
        radius, speed, wind, gamma = 2e10, 1e6, 4e5, 5e-12
        target, horizon = replay.AU_M, 3e5

        def rhs(_, state):
            relative_speed = state[1] - wind
            return [state[1], -gamma * relative_speed * abs(relative_speed)]

        def crossing(_, state):
            return state[0] - target

        # Integrate the governing ODE independently of the core analytic drag
        # expression. Verify both the crossing time and the speed at the root.
        reference = solve_ivp(rhs, (0, horizon), [radius, speed], events=crossing,
                              rtol=1e-11, atol=1e-7, max_step=600)
        reply = replay.Core(CORE, []).batch([['A', radius, speed, wind, gamma,
                                            target, horizon]])[0]
        self.assertTrue(reference.success)
        self.assertEqual(reply[0], 'ARRIVAL')
        self.assertAlmostEqual(float(reply[1]), reference.t_events[0][0], places=3)
        self.assertAlmostEqual(float(reply[2]), reference.y_events[0][0][1], places=2)

    def test_noisy_two_views_and_causal_cut(self):
        event = synthetic(noise=True)
        result = replay.replay(event, replay.utc(timestamp(11400)), {'C'}, CORE, Path.cwd())
        any_ignored = result['fit']['excluded']
        self.assertEqual(any_ignored['not_yet_available'], 1)
        self.assertGreater(any_ignored['withheld_view'], 0)
        np.testing.assert_allclose(result['fit']['parameters'], TRUTH, atol=0.005, rtol=0)
        self.assertGreater(result['fit']['reduced_chi_squared'], 0.5)
        self.assertLess(result['fit']['reduced_chi_squared'], 1.5)

    def test_production_cli_maps_and_one_au_crossing(self):
        event = synthetic()
        with tempfile.TemporaryDirectory() as temp:
            directory = Path(temp)
            data, output = directory / 'event.json', directory / 'out.json'
            data.write_text(json.dumps(event))
            subprocess.run([sys.executable, str(Path(replay.__file__)), str(data),
                            '--core', str(CORE), '--issue', timestamp(11400),
                            '--output', str(output)], check=True)
            report = json.loads(output.read_bytes())
        apex, outside = report['map']
        expected_height = 8 + 7.5 * 3 + 0.125 * 9
        expected_v = 7.5 + 0.25 * 3
        np.testing.assert_allclose(report['fit']['parameters'], TRUTH, atol=1e-8, rtol=0)
        self.assertAlmostEqual(np.linalg.norm(apex['nominal_geometry']['position_hci_m'])
                               / replay.SOLAR_RADIUS_M, expected_height, places=8)
        self.assertAlmostEqual(apex['nominal_geometry']['normal_speed_m_s']
                               * 3600 / replay.SOLAR_RADIUS_M, expected_v, places=8)
        self.assertEqual(apex['statuses'], {replay.SHOCK: 16})
        self.assertEqual(outside['statuses'], {'outside-front-support': 16})
        self.assertGreater(apex['compression_given_solved_p05_p50_p95'][1], 1)
        self.assertLess(apex['compression_given_solved_p05_p50_p95'][1], 4)
        # Every usable acquisition epoch gets the same physical-ray map. The
        # independently manufactured radius and speed must evolve at each one.
        self.assertEqual([frame['elapsed_s'] for frame in report['evolution']],
                         [0, 1188, 4320, 7200, 10800])
        self.assertEqual(report['evolution'][-1]['map'], report['map'])
        for frame in report['evolution']:
            hour = frame['elapsed_s'] / 3600
            self.assertEqual(frame['source'], 'fitted-history')
            self.assertEqual(frame['map'][0]['ray_hci'], apex['ray_hci'])
            self.assertAlmostEqual(np.linalg.norm(frame['map'][0]['nominal_geometry']['position_hci_m'])
                                   / replay.SOLAR_RADIUS_M, 8 + 7.5 * hour + 0.125 * hour**2,
                                   places=7)
            self.assertAlmostEqual(frame['map'][0]['nominal_geometry']['normal_speed_m_s']
                                   * 3600 / replay.SOLAR_RADIUS_M, 7.5 + 0.25 * hour, places=7)
        self.assertNotEqual(report['evolution'][0]['map'][0]['fast_mach_given_solved_p05_p50_p95'],
                            apex['fast_mach_given_solved_p05_p50_p95'])
        # With zero drag, the independent 1-D nose crossing is analytic.
        # The physical classifier at the crossing still uses the PARKER state.
        eta = (3 + (replay.AU_M / replay.SOLAR_RADIUS_M - expected_height)
               / expected_v) * 3600
        arrival = report['observers']['nose']
        self.assertEqual(arrival['fast_shock_fraction_bounds'], [1, 1])
        self.assertLess(abs(arrival['arrival_elapsed_s_given_solved_p05_p50_p95'][1] - eta), 60)
        self.assertEqual(report['observers']['far_side']['fast_shock_fraction_bounds'], [0, 0])

    def test_requested_evolution_matches_drag_ode_without_observers(self):
        event = synthetic()
        event.pop('observers_hci_m')
        last_s, late_s = 10800, 18000
        event['evolution_elapsed_s'] = [late_s, 1800.5, last_s, last_s]
        wind, gamma = 4e5, 5e-12
        event['propagation']['drag_m_inv'] = gamma
        result = replay.replay(event, replay.utc(timestamp(11400)), set(), CORE, Path.cwd())
        early, handoff, outer = result['evolution']
        self.assertEqual([f['elapsed_s'] for f in result['evolution']], [1800.5, last_s, late_s])
        self.assertEqual(result['observers'], {})
        self.assertEqual(early['source'], 'fitted-history')
        self.assertEqual(outer['source'], 'quadratic-drag')

        def rhs(_, state):
            relative_speed = state[1] - wind
            return [state[1], -gamma * relative_speed * abs(relative_speed)]

        # The independent fit truth fixes the drag initial condition. Solve
        # its ODE at the requested time to check both the handoff and evolution.
        initial = [(8 + 7.5 * 3 + 0.125 * 9) * replay.SOLAR_RADIUS_M,
                   (7.5 + 0.25 * 3) * replay.SOLAR_RADIUS_M / 3600]
        solution = solve_ivp(rhs, (0, late_s - last_s), initial,
                             rtol=1e-11, atol=1e-7, max_step=600)
        self.assertTrue(solution.success)
        np.testing.assert_allclose([handoff['apex_radius_m'], handoff['apex_speed_m_s']],
                                   initial, rtol=1e-9)
        np.testing.assert_allclose([outer['apex_radius_m'],
                                    outer['map'][0]['nominal_geometry']['normal_speed_m_s']],
                                   solution.y[:, -1], rtol=1e-9)
        self.assertEqual(outer['map'][0]['statuses'], {replay.SHOCK: 16})

    def test_evolution_rejects_unqualified_extrapolation(self):
        event = synthetic()
        event.pop('observers_hci_m')
        event['evolution_elapsed_s'] = [10801]
        propagation = event.pop('propagation')
        issue = replay.utc(timestamp(11400))
        with self.assertRaisesRegex(ValueError, 'requires configured propagation'):
            replay.replay(event, issue, set(), CORE, Path.cwd())
        event['propagation'] = propagation
        event['evolution_elapsed_s'] = [200000]
        with self.assertRaisesRegex(ValueError, 'exceeds the declared drag horizon'):
            replay.replay(event, issue, set(), CORE, Path.cwd())
        event['evolution_elapsed_s'] = [10801]
        event['propagation']['minimum_handoff_apex_m'] = 40 * replay.SOLAR_RADIUS_M
        result = replay.replay(event, issue, set(), CORE, Path.cwd())
        future = result['evolution'][0]
        self.assertEqual(future['trajectory_status'], 'BELOW_DRAG_HANDOFF')
        self.assertIsNone(future['apex_radius_m'])
        self.assertEqual(future['map'][0]['statuses'], {'BELOW_DRAG_HANDOFF': 16})
        self.assertEqual(future['map'][0]['unresolved_fraction'], 1)

    def test_missing_ambient_and_censored_arrival(self):
        event = synthetic()
        event['ambient_members'] = []
        result = replay.replay(event, replay.utc(timestamp(11400)), set(), CORE, Path.cwd())
        self.assertEqual(result['map'][0]['unresolved_fraction'], 1)
        self.assertIsNotNone(result['map'][0]['nominal_geometry']['position_hci_m'])
        observer = result['observers']['nose']
        self.assertIsNotNone(observer['front_crossing_elapsed_s_p05_p50_p95'])
        self.assertIsNone(observer['arrival_elapsed_s_given_solved_p05_p50_p95'])
        event['propagation']['horizon_s'] = 100
        result = replay.replay(event, replay.utc(timestamp(11400)), set(), CORE, Path.cwd())
        observer = result['observers']['nose']
        self.assertEqual(observer['statuses'], {'BEYOND_FORECAST_HORIZON': 16})
        self.assertIsNone(observer['front_crossing_elapsed_s_p05_p50_p95'])

    def test_upstream_ensemble_keeps_missing_member_weight(self):
        with tempfile.TemporaryDirectory() as temp:
            directory = Path(temp)
            shifted = directory / 'future.event'
            shutil.copytree(AMBIENT.parent / 'assets', directory / 'assets')
            # This manufactured member is already available at issue and its
            # validity begins at +24 h. At the map epoch it has no coverage.
            # The observer crossing is later, inside both model time ranges.
            shifted.write_text(AMBIENT.read_text().replace(REFERENCE, '2026-10-07T00:00:00Z'))
            event = synthetic()
            event['ambient_members'].extend([
                {'event_file': shifted.name, 'available_utc': REFERENCE},
                {'event_file': 'not-yet-received.event', 'available_utc': timestamp(11500)},
            ])
            result = replay.replay(event, replay.utc(timestamp(11400)), set(), CORE, directory)
        self.assertEqual(result['ambient_members'][1]['event_to_ambient_offset_s'], -86400)
        self.assertEqual(len(result['ambient_not_available_at_issue']), 1)
        report = result['map'][0]
        self.assertEqual(report['statuses'], {replay.SHOCK: 16, 'ambient-unavailable': 16})
        self.assertEqual(report['fast_shock_fraction_bounds'], [0.5, 1])
        self.assertEqual(report['unresolved_fraction'], 0.5)
        self.assertEqual(result['observers']['nose']['statuses'], {replay.SHOCK: 32})

    def test_rank_deficiency_and_atomic_output(self):
        event = synthetic()
        distinct = {}
        for track in event['tracks']:
            if track['view'] == 'A':
                distinct.setdefault(track['elapsed_s'], track)
        event['tracks'] = [track for track in distinct.values() for _ in range(8)]
        with self.assertRaisesRegex(ValueError, 'underdetermined|no regular'):
            replay.fit_tracks(event, replay.utc(timestamp(11400)), set())
        with tempfile.TemporaryDirectory() as temp:
            directory = Path(temp)
            data, output = directory / 'event.json', directory / 'out.json'
            data.write_text(json.dumps(event))
            output.write_text('previous receipt')
            rejected = subprocess.run([sys.executable, str(Path(replay.__file__)), str(data),
                                       '--core', str(CORE), '--issue', timestamp(11400),
                                       '--output', str(output)], capture_output=True, text=True)
            self.assertEqual(rejected.returncode, 2)
            self.assertEqual(output.read_text(), 'previous receipt')


if __name__ == '__main__':
    unittest.main(verbosity=2)
