#!/usr/bin/env python3
"""Pure tests for Stage-3 particle-time adaptive radial bins."""

from __future__ import annotations

import math
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "srcHartley2" / "visualization"))

from adaptive_radial_bins import (  # noqa: E402
    Cell,
    PlotError,
    load_cells,
    main,
    merge_cells,
    merge_replicate_cells,
)


class AdaptiveBinTest(unittest.TestCase):
    def test_threshold_support_and_integrated_number(self) -> None:
        cells = [
            Cell(i, i + 1.5, i + 1.0, i + 2.0, float(i + 1),
                 density_m3=(-1.0 if i == 2 else 10.0 * (i + 1)),
                 occupancy=value)
            for i, value in enumerate((0.0, 1.0, 2.0, 0.0, 0.5))
        ]
        bins = merge_cells(cells, sampling_length=10, minimum_exposure=20.0)
        self.assertEqual([(row["first_cell"], row["last_cell"]) for row in bins],
                         [(0, 2), (3, 4)])
        self.assertEqual(bins[0]["particle_time_exposure"], 30.0)
        self.assertEqual(bins[1]["status"], "final_underexposed")
        self.assertEqual(bins[0]["nonpositive_cell_count"], 1)
        source_integrated = sum(cell.density_m3 * cell.volume_m3 for cell in cells)
        merged_integrated = sum(row["integrated_number"] for row in bins)
        self.assertTrue(math.isclose(source_integrated, merged_integrated))
        self.assertEqual(bins[0]["r_outer_m"], bins[1]["r_inner_m"])

    def test_all_zero_exposure_is_retained(self) -> None:
        cells = [Cell(i, i + 1.5, i + 1.0, i + 2.0, 1.0, 0.0, 0.0)
                 for i in range(3)]
        bins = merge_cells(cells, 100, 50.0)
        self.assertEqual(len(bins), 1)
        self.assertEqual(bins[0]["status"], "zero_exposure")
        self.assertEqual(bins[0]["cell_count"], 3)

    def test_negative_exposure_is_fatal(self) -> None:
        cells = [Cell(0, 1.5, 1.0, 2.0, 1.0, 1.0, -1.0)]
        with self.assertRaises(ValueError):
            merge_cells(cells, 10, 1.0)

    def test_independent_seed_uncertainty_uses_common_support(self) -> None:
        replicas = []
        for seed_density in (8.0, 10.0, 12.0, 14.0):
            replicas.append(
                [
                    Cell(
                        i,
                        i + 1.5,
                        i + 1.0,
                        i + 2.0,
                        1.0,
                        seed_density + i,
                        occupancy=(0.5 if i == 0 else 1.5),
                    )
                    for i in range(2)
                ]
            )
        bins = merge_replicate_cells(replicas, 10, 20.0)
        self.assertEqual(len(bins), 1)
        self.assertEqual(bins[0]["first_cell"], 0)
        self.assertEqual(bins[0]["last_cell"], 1)
        self.assertEqual(bins[0]["uncertainty_sample_count"], 4)
        self.assertEqual(bins[0]["decision_sample_status"], "ready")
        self.assertTrue(
            math.isclose(bins[0]["number_density_m-3"], 11.5)
        )
        self.assertIsNotNone(
            bins[0]["density_sample_standard_deviation_m-3"]
        )

    def test_replicate_mesh_mismatch_is_fatal(self) -> None:
        one = [Cell(0, 1.5, 1.0, 2.0, 1.0, 1.0, 1.0)]
        other = [Cell(0, 1.6, 1.0, 2.0, 1.0, 1.0, 1.0)]
        with self.assertRaises(ValueError):
            merge_replicate_cells([one, other], 10, 1.0)

    def test_center_cell_csv_joins_by_exact_index(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            measures = root / "cell_measures.csv"
            samples = root / "stage3_cells.C2.out=3.csv"
            measures.write_text(
                "cell_index,radius_m,r_inner_m,r_outer_m,measure_m3\n"
                "0,1.5,1,2,10\n1,2.5,2,3,20\n",
                encoding="utf-8",
            )
            samples.write_text(
                "cell_index,radius_m,r_inner_m,r_outer_m,measure_m3,number_density_m-3,mean_model_particles,sampling_length,output_index,random_seed,random_seed_configured,residual_sampling,mpi_ranks\n"
                "0,1.5,1,2,10,4,0.5,100,3,101,1,SYSTEMATIC,4\n"
                "1,2.5,2,3,20,8,2.5,100,3,101,1,SYSTEMATIC,4\n",
                encoding="utf-8",
            )
            cells = load_cells(samples, measures)
            self.assertEqual(len(cells), 2)
            self.assertEqual(cells[0].density_m3, 4.0)
            self.assertEqual(cells[1].occupancy, 2.5)

    def test_optional_block_policy_metadata_is_preserved(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            measures = root / "cell_measures.csv"
            samples = root / "stage3_cells.C2.out=3.csv"
            measures.write_text(
                "cell_index,radius_m,r_inner_m,r_outer_m,measure_m3\n"
                "0,1.5,1,2,10\n1,2.5,2,3,20\n",
                encoding="utf-8",
            )
            samples.write_text(
                "cell_index,block_index,local_cell_index,block_cell_count,"
                "radius_m,r_inner_m,r_outer_m,measure_m3,number_density_m-3,"
                "mean_model_particles,population_target_min_per_cell,"
                "population_target_max_per_cell,population_locked,"
                "population_class,population_policy_source,sampling_length,"
                "output_index,random_seed,random_seed_configured,"
                "residual_sampling,mpi_ranks\n"
                "0,0,0,2,1.5,1,2,10,4,0.5,50,200,0,2,3,100,3,101,1,SYSTEMATIC,4\n"
                "1,0,1,2,2.5,2,3,20,8,2.5,50,200,0,2,3,100,3,101,1,SYSTEMATIC,4\n",
                encoding="utf-8",
            )
            cells = load_cells(samples, measures)
            self.assertEqual(cells[0].block_index, 0)
            self.assertEqual(cells[1].local_cell_index, 1)
            self.assertEqual(cells[1].block_cell_count, 2)
            self.assertEqual(cells[0].population_target_min, 50)
            self.assertEqual(cells[0].population_target_max, 200)
            self.assertEqual(cells[0].population_class, 2)
            self.assertEqual(cells[0].population_policy_source, 3)

    def test_block_policy_indices_must_match_cell_order(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            measures = root / "cell_measures.csv"
            samples = root / "stage3_cells.C2.out=3.csv"
            measures.write_text(
                "cell_index,radius_m,r_inner_m,r_outer_m,measure_m3\n"
                "0,1.5,1,2,10\n",
                encoding="utf-8",
            )
            samples.write_text(
                "cell_index,block_index,local_cell_index,block_cell_count,"
                "radius_m,r_inner_m,r_outer_m,measure_m3,number_density_m-3,"
                "mean_model_particles,population_target_min_per_cell,"
                "population_target_max_per_cell,population_locked,"
                "population_class,population_policy_source,sampling_length,"
                "output_index,random_seed,random_seed_configured,"
                "residual_sampling,mpi_ranks\n"
                "0,1,0,2,1.5,1,2,10,4,0.5,50,200,0,2,3,100,3,101,1,SYSTEMATIC,4\n",
                encoding="utf-8",
            )
            with self.assertRaises(PlotError):
                load_cells(samples, measures)

    def test_command_line_path_writes_both_output_formats(self) -> None:
        """Exercise the complete CLI path so reporting-only failures are caught."""
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            measures = root / "cell_measures.csv"
            samples = root / "stage3_cells.C2.out=3.csv"
            output = root / "adaptive"
            measures.write_text(
                "cell_index,radius_m,r_inner_m,r_outer_m,measure_m3\n"
                "0,1.5,1,2,10\n1,2.5,2,3,20\n",
                encoding="utf-8",
            )
            samples.write_text(
                "cell_index,radius_m,r_inner_m,r_outer_m,measure_m3,number_density_m-3,mean_model_particles,sampling_length,output_index,random_seed,random_seed_configured,residual_sampling,mpi_ranks\n"
                "0,1.5,1,2,10,4,0,100,3,101,1,SYSTEMATIC,4\n"
                "1,2.5,2,3,20,8,0,100,3,101,1,SYSTEMATIC,4\n",
                encoding="utf-8",
            )
            self.assertEqual(
                main(
                    [
                        str(samples),
                        "--cell-measures",
                        str(measures),
                        "--sampling-length",
                        "100",
                        "--minimum-exposure",
                        "1",
                        "-o",
                        str(output),
                    ]
                ),
                0,
            )
            self.assertTrue(output.with_suffix(".csv").is_file())
            self.assertTrue(output.with_suffix(".json").is_file())

    def test_cli_rejects_duplicate_seed_replicates(self) -> None:
        """Four path arguments are not four independent statistical samples."""
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            measures = root / "cell_measures.csv"
            measures.write_text(
                "cell_index,radius_m,r_inner_m,r_outer_m,measure_m3\n"
                "0,1.5,1,2,10\n",
                encoding="utf-8",
            )
            tables = []
            for suffix in ("a", "b"):
                run = root / suffix
                run.mkdir()
                table = run / "stage3_cells.C2.out=3.csv"
                table.write_text(
                    "cell_index,radius_m,r_inner_m,r_outer_m,measure_m3,number_density_m-3,mean_model_particles,sampling_length,output_index,random_seed,random_seed_configured,residual_sampling,mpi_ranks\n"
                    "0,1.5,1,2,10,4,1,100,3,101,1,SYSTEMATIC,4\n",
                    encoding="utf-8",
                )
                tables.append(table)
            with self.assertRaises(PlotError):
                main(
                    [
                        str(tables[0]),
                        "--replicate",
                        str(tables[1]),
                        "--cell-measures",
                        str(measures),
                        "--sampling-length",
                        "100",
                        "--minimum-exposure",
                        "1",
                        "-o",
                        str(root / "adaptive"),
                    ]
                )

    def test_cli_rejects_sampling_length_override_mismatch(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            measures = root / "cell_measures.csv"
            samples = root / "stage3_cells.C2.out=3.csv"
            measures.write_text(
                "cell_index,radius_m,r_inner_m,r_outer_m,measure_m3\n"
                "0,1.5,1,2,10\n",
                encoding="utf-8",
            )
            samples.write_text(
                "cell_index,radius_m,r_inner_m,r_outer_m,measure_m3,number_density_m-3,mean_model_particles,sampling_length,output_index,random_seed,random_seed_configured,residual_sampling,mpi_ranks\n"
                "0,1.5,1,2,10,4,1,100,3,101,1,SYSTEMATIC,4\n",
                encoding="utf-8",
            )
            with self.assertRaises(PlotError):
                main(
                    [
                        str(samples),
                        "--cell-measures",
                        str(measures),
                        "--sampling-length",
                        "50",
                        "--minimum-exposure",
                        "1",
                        "-o",
                        str(root / "adaptive"),
                    ]
                )


if __name__ == "__main__":
    unittest.main(verbosity=2)
