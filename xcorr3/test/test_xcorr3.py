"""Portable launcher contract tests; run in PowerShell or Linux with Python 3.

Backend selection is mocked; subprocess tests execute Python fixtures, not SAR
processing. These tests do not claim Linux/CUDA numerical validation.
"""
import contextlib
import importlib.machinery
import importlib.util
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import Mock, patch

ENTRY = Path(__file__).resolve().parents[1] / 'xcorr3'
loader = importlib.machinery.SourceFileLoader('xcorr3_api', str(ENTRY))
spec = importlib.util.spec_from_loader(loader.name, loader)
api = importlib.util.module_from_spec(spec)
loader.exec_module(api)


class LauncherTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix='xcorr3-tests-')
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.config = self.root / 'config.txt'
        self.config.write_text('xcorr_backend = mt\nxcorr_nproc = 3\n', encoding='utf-8')

    def invoke(self, args, resolve=True, result=0):
        out, err = io.StringIO(), io.StringIO()
        with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err), \
             patch.object(api, 'resolve_backend', return_value='/fake/backend') as finder, \
             patch.object(api.subprocess, 'run', return_value=subprocess.CompletedProcess([], result)) as run:
            if not resolve:
                finder.side_effect = ValueError('backend executable not found')
            try:
                status = api.main(args)
            except SystemExit as error:
                status = error.code
        return status, out.getvalue(), err.getvalue(), run, finder

    def dry(self, args):
        status, out, err, run, _ = self.invoke(['--dry-run', *args])
        self.assertEqual(status, 0, err)
        run.assert_not_called()
        return json.loads(out)

    def test_missing_input_prints_complete_usage(self):
        for args in ([], ['master.PRM'], ['--backend', 'cc'], ['--config', str(self.config)]):
            with self.subTest(args=args):
                status, _, err, run, finder = self.invoke(args)
                self.assertEqual(status, 2)
                self.assertIn('Required inputs:', err)
                self.assertIn('Examples:', err)
                run.assert_not_called(); finder.assert_not_called()

    def test_help_without_installation(self):
        result = subprocess.run([sys.executable, str(ENTRY), '--help'], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0)
        self.assertIn('master.PRM secondary.PRM', result.stdout)

    def test_real_cli_missing_arguments(self):
        result = subprocess.run([sys.executable, str(ENTRY), 'only.PRM'], capture_output=True, text=True)
        self.assertEqual(result.returncode, 2)
        self.assertIn('Examples:', result.stderr)

    def test_default_original(self):
        data = self.dry(['a.PRM', 'b.PRM'])
        self.assertEqual(data['backend'], 'original')
        self.assertEqual(data['arguments'], ['a.PRM', 'b.PRM'])

    def test_backend_before_and_after_inputs(self):
        for name in api.BACKENDS:
            for args in (['--backend', name, 'a', 'b'], ['a', 'b', '--backend', name]):
                self.assertEqual(self.dry(args)['backend'], name)

    def test_config_selection(self):
        data = self.dry(['--config', str(self.config), 'a', 'b'])
        self.assertEqual(data['backend'], 'mt')
        self.assertEqual(data['workers'], '3')
        self.assertEqual(data['arguments'], ['a', 'b', '-nproc', '3'])

    def test_cli_backend_overrides_config(self):
        data = self.dry(['--config', str(self.config), '--backend', 'cc', 'a', 'b'])
        self.assertEqual(data['backend'], 'cc')
        self.assertNotIn('-nproc', data['arguments'])

    def test_worker_override_and_auto(self):
        for flag in ('-nproc', '--nproc'):
            for count in ('2', 'auto'):
                args = self.dry(['--config', str(self.config), 'a', 'b', flag, count])['arguments']
                self.assertEqual(args, ['a', 'b'] + (['-nproc', '2'] if count == '2' else []))

    def test_original_proposed_separator_syntax(self):
        args = self.dry(['--backend', 'mt', '--nproc', '6', '--', 'xcorr', 'a', 'b'])['arguments']
        self.assertEqual(args, ['a', 'b', '-nproc', '6'])

    def test_native_worker_after_separator(self):
        data = self.dry(['--backend', 'mt', '--', 'xcorr', 'a', 'b', '-nproc', '4'])
        self.assertEqual(data['arguments'], ['a', 'b', '-nproc', '4'])

    def test_invalid_workers(self):
        for value in ('0', '-1', 'abc', '1.2', '2147483648', '', '６'):
            status, _, _, run, _ = self.invoke(['--backend', 'mt', 'a', 'b', '-nproc', value])
            self.assertEqual(status, 2); run.assert_not_called()

    def test_missing_option_value(self):
        for option in ('--backend', '--config', '--nproc', '-nproc'):
            status, _, err, run, _ = self.invoke(['a', 'b', option])
            self.assertEqual(status, 2); self.assertIn('Examples:', err); run.assert_not_called()

    def test_duplicate_workers(self):
        for extra in (['-nproc', '2', '-nproc', '3'], ['--nproc', '2', '-nproc', '3']):
            self.assertEqual(self.invoke(['--backend', 'mt', 'a', 'b', *extra])[0], 2)

    def test_workers_rejected_for_other_backends(self):
        for backend in ('original', 'cc'):
            for flag in ('-nproc', '--nproc'):
                self.assertEqual(self.invoke(['--backend', backend, 'a', 'b', flag, '2'])[0], 2)

    def test_invalid_and_duplicate_config(self):
        for text in ('xcorr_backend=other', 'xcorr_nproc=0', 'xcorr_backend mt',
                     'xcorr_backend extra=mt', 'xcorr_backend=mt\nxcorr_backend=cc'):
            self.config.write_text(text, encoding='utf-8')
            self.assertEqual(self.invoke(['--config', str(self.config), 'a', 'b'])[0], 2)

    def test_config_comments_bom_and_unrelated_fields(self):
        self.config.write_text('\ufeffproc_stage=1\nxcorr_backend=mt # comment\nxcorr_nproc=auto\n', encoding='utf-8')
        self.assertEqual(self.dry(['--config', str(self.config), 'a', 'b'])['arguments'], ['a', 'b'])

    def test_empty_config_defaults(self):
        self.config.write_text('xcorr_backend=\nxcorr_nproc=\n', encoding='utf-8')
        self.assertEqual(self.dry(['--config', str(self.config), 'a', 'b'])['backend'], 'original')

    def test_missing_config(self):
        self.assertEqual(self.invoke(['--config', str(self.root/'absent'), 'a', 'b'])[0], 2)

    def test_cc_unsupported_modes(self):
        for mode in ('-time', '-time4', '-real', '-grid', '-af'):
            self.assertEqual(self.invoke(['--backend', 'cc', 'a', 'b', mode])[0], 2)

    def test_cc_only_options(self):
        for backend in ('original', 'mt'):
            for flag in ('-geocode', '-psnr', '-nyquist_split'):
                self.assertEqual(self.invoke(['--backend', backend, 'a', 'b', flag])[0], 2)
        self.assertIn('-geocode', self.dry(['--backend', 'cc', 'a', 'b', '-geocode'])['arguments'])

    def test_explicit_frequency_mode(self):
        for backend in api.BACKENDS:
            args=self.dry(['--backend', backend, 'a', 'b', '-freq', '-nx', '20'])['arguments']
            self.assertEqual(args, ['a','b']+([] if backend=='cc' else ['-freq'])+['-nx','20'])

    def test_parameter_passthrough_and_space_paths(self):
        args = ['a = image.PRM', 'b image.PRM', '-nx', '20', '-ny', '50', '-interp', '3', '-nointerp']
        self.assertEqual(self.dry(['--backend', 'cc', *args])['arguments'], args)

    def test_no_fallback_on_missing_executable(self):
        status, _, err, run, finder = self.invoke(['--backend', 'cc', 'a', 'b'], resolve=False)
        self.assertEqual(status, 2); run.assert_not_called(); finder.assert_called_once_with('cc')
        self.assertIn('not found', err)

    def test_child_status_and_signal(self):
        for status, expected in ((0, 0), (7, 7), (-15, 143)):
            result = self.invoke(['--backend', 'cc', 'a', 'b'], result=status)
            self.assertEqual(result[0], expected); self.assertEqual(result[3].call_count, 1)

    def test_actual_process_and_argument_boundaries(self):
        script = self.root / 'master = fixture.PRM'
        script.write_text('import json,sys\nprint(json.dumps(sys.argv[1:]))\nsys.exit(7)\n', encoding='utf-8')
        driver = self.root / 'driver.py'
        driver.write_text('import runpy,sys\na=runpy.run_path(sys.argv[1])\n'
            'sys.exit(a["dispatch"](sys.argv[2:],"cc","auto",sys.executable))\n', encoding='utf-8')
        args = ['secondary ; name.PRM', '-nx', '20']
        result = subprocess.run([sys.executable, str(driver), str(ENTRY), str(script), *args], capture_output=True, text=True)
        self.assertEqual(result.returncode, 7)
        self.assertEqual(json.loads(result.stdout), args)

    def test_backend_discovery(self):
        # Isolate discovery from binaries already built in the source tree.
        with patch.object(api, '__file__', str(self.root/'xcorr3')), patch.object(api.shutil, 'which', return_value=sys.executable):
            for backend in api.BACKENDS:
                self.assertEqual(api.resolve_backend(backend), str(Path(sys.executable).resolve()))
        with patch.object(api, '__file__', str(self.root/'xcorr3')), patch.object(api.shutil, 'which', return_value=None):
            with self.assertRaises(ValueError): api.resolve_backend('cc')

    def test_installed_sibling_precedes_path(self):
        binary = self.root/'xcorr_cc'
        binary.write_text('fixture', encoding='utf-8'); binary.chmod(0o755)
        with patch.object(api, '__file__', str(self.root/'xcorr3')), patch.object(api.shutil, 'which', return_value=None):
            self.assertEqual(api.resolve_backend('cc'), str(binary.resolve()))

    def test_environment_unchanged(self):
        old = dict(os.environ)
        self.invoke(['--backend', 'mt', 'a', 'b'])
        self.assertEqual(dict(os.environ), old)

    def test_workflow_contract(self):
        with patch.object(api.shutil, 'which', return_value='/fake/csh'):
            for command in (['p2p_processing.csh', 'S1_TOPS'], ['csh', 'align_tops.csh']):
                with self.assertRaises(ValueError): api.validate_workflow(command, 'cc')
            with self.assertRaises(ValueError): api.validate_workflow(['/usr/bin/xcorr'], 'mt')
            with patch.dict(os.environ, XCORR3_ACTIVE='1'):
                with self.assertRaises(ValueError): api.validate_workflow(['p2p_processing.csh'], 'mt')

    def test_workflow_dry_run(self):
        with patch.object(api.shutil, 'which', return_value='/fake/workflow'):
            data = self.dry(['--backend', 'mt', '--nproc', '2', '--run', 'workflow', '-flag'])
        self.assertEqual(data['mode'], 'workflow'); self.assertEqual(data['arguments'], ['workflow', '-flag'])
        self.assertEqual(data['workers'], '2')

    def test_workflow_missing_or_mixed_inputs(self):
        for args in (['--run'], ['a', 'b', '--run', 'workflow']):
            self.assertEqual(self.invoke(args)[0], 2)

    def test_workflow_shim_dispatch_and_cleanup(self):
        # Exercise the generated shim on Windows via Python, not a csh claim.
        script = self.root / 'input with spaces.PRM'
        script.write_text('import sys\nsys.exit(9 if sys.argv[1:]==["b.PRM","-nproc","2"] else 1)\n', encoding='utf-8')
        actual_run = subprocess.run
        observed = {}
        def launch(command, env):
            directory = Path(env['PATH'].split(os.pathsep)[0])
            observed['directory'] = directory
            self.assertEqual(command, ['workflow'])
            self.assertEqual(env['XCORR3_ACTIVE'], '1')
            self.assertEqual(env['XCORR3_EXECUTABLE'], sys.executable)
            child = actual_run([sys.executable, str(directory/'xcorr'), str(script), 'b.PRM', '-nproc', '2'],
                               env=env, capture_output=True, text=True)
            self.assertEqual(child.returncode, 9)
            # Simulate a shell swallowing the backend error with a final true.
            return subprocess.CompletedProcess(command, 0)
        proxy = Mock(wraps=os); proxy.name = 'posix'; proxy.pathsep = os.pathsep
        with patch.object(api, 'os', proxy), patch.object(api.shutil, 'which', return_value='/fake/workflow'), \
             patch.object(api.subprocess, 'run', side_effect=launch):
            status = api.run_workflow(['workflow'], 'mt', '6', sys.executable)
        self.assertEqual(status, 9)
        self.assertFalse(observed['directory'].exists())

    @unittest.skipUnless(os.name == 'nt', 'Windows-specific platform guard')
    def test_windows_workflow_rejected(self):
        with self.assertRaisesRegex(ValueError, 'POSIX'):
            api.run_workflow(['workflow'], 'cc', 'auto', '/fake/backend')


if __name__ == '__main__':
    unittest.main(verbosity=2)
