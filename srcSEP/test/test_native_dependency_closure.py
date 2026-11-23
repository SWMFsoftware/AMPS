#!/usr/bin/env python3
"""Ensure stale copied PIC templates cannot qualify native mover tests."""
import os
from pathlib import Path
import sys
import subprocess
import tempfile
import unittest

sys.dont_write_bytecode = True
from run_native_standalone_integration import pic_dependency_closure


class DependencyControls(unittest.TestCase):
    def assert_executable_dependencies(self, text):
        # Execute the REAL Make dependency stanza in a small owned fixture.
        # Header changes must rebuild both main and debug archive members;
        # merely recognizing stale objects in a later auditor is insufficient.
        begin=min(text.index('MAINOBJ ='),text.index('# Bootstrap dependency files'))
        end=text.index('\nlib:',begin)
        stanza=text[begin:end]
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);(root/'debug').mkdir()
            for name in ('main.o','debug/sep_debug_fieldline_datum.o','lib.o','picGlobal.dfn'):
                (root/name).write_text('owned fixture')
                os.utime(root/name,ns=(2000000000,2000000000))
            os.utime(root/'picGlobal.dfn',ns=(3000000000,3000000000))
            for name in ('main','debug/sep_debug_fieldline_datum'):
                (root/(name+'.d')).write_text(name+'.o: picGlobal.dfn\n')
            (root/'lib.d').write_text('lib.o:\n')
            (root/'fixture.mk').write_text('MAINLIBOBJ = lib.o\n'+stanza+
                '\nmain.o:\n\t@echo REBUILD_MAIN\ndebug/sep_debug_fieldline_datum.o:\n\t@echo REBUILD_DEBUG\n')
            result=subprocess.run(['make','-n','-f','fixture.mk','main.o','debug/sep_debug_fieldline_datum.o'],
                                  cwd=root,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
            self.assertEqual(result.returncode,0,result.stdout)
            self.assertIn('REBUILD_MAIN',result.stdout)
            self.assertIn('REBUILD_DEBUG',result.stdout)

    def test_main_and_debug_follow_real_header_dependencies(self):
        self.assert_executable_dependencies((Path(__file__).resolve().parents[1]/'makefile').read_text())

    def fixture(self, root):
        pic, main = root/"build/pic", root/"build/main"
        pic.mkdir(parents=True)
        (main/"util").mkdir(parents=True)
        header, source, obj = pic/"pic.h", main/"util/probe.cpp", main/"util/probe.o"
        for path in (header, source, obj):
            path.write_text("test dependency")
        (main/"util/probe.d").write_text(
            "util/probe.o: util/probe.cpp " + str(header) + "\n\n" + str(header) + ":\n")
        for path in (header, source):
            os.utime(path, ns=(1000000000, 1000000000))
        os.utime(obj, ns=(2000000000, 2000000000))
        return header, source, obj

    def test_fresh_nested_compilation_directory(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            header, source, obj = self.fixture(root)
            paths = pic_dependency_closure(root)
            self.assertIn(header, paths)
            self.assertIn(source, paths)
            self.assertIn(obj, paths)

    def test_newer_pic_header_rejects_old_template_object(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            header, _, _ = self.fixture(root)
            os.utime(header, ns=(3000000000, 3000000000))
            with self.assertRaisesRegex(ValueError, "stale native PIC"):
                pic_dependency_closure(root)

    def test_other_build_header_cannot_supply_pic_authority(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            header, _, _ = self.fixture(root)
            other = root/"old-build/pic.h"
            other.parent.mkdir()
            other.write_bytes(header.read_bytes())
            os.utime(other, ns=(1000000000, 1000000000))
            depfile = root/"build/main/util/probe.d"
            depfile.write_text(depfile.read_text().replace(str(header), str(other)))
            with self.assertRaisesRegex(ValueError, "another configured build"):
                pic_dependency_closure(root)


if __name__ == "__main__":
    unittest.main()
