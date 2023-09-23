import unittest
from unittest.mock import patch, MagicMock
import sys, os

# 假设你的arg_parser在scc包里面
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from scc.arg_parser import get_args_and_check_file

# 模拟 argparse 的 Namespace 对象
class MockArgs:
    def __init__(self, **entries):
        self.__dict__.update(entries)

class TestArgParser(unittest.TestCase):

    @patch('os.path.exists')
    @patch('argparse.ArgumentParser.parse_args')
    def test_file_existence(self, mock_parse_args, mock_exists):
        mock_parse_args.return_value = MockArgs(
            bam1='bam1_file',
            bam2='bam2_file',
            ref1='ref1_file',
            ref2='ref2_file',
            snp='snp_file',
            outprefix='outprefix'
        )

        # 模拟文件存在
        mock_exists.return_value = True

        args = get_args_and_check_file()

        # 断言文件存在检查被调用了
        mock_exists.assert_called()

        # 断言返回的args是正确的
        self.assertEqual(args.bam1, 'bam1_file')

    @patch('os.path.exists')
    @patch('argparse.ArgumentParser.parse_args')
    def test_file_nonexistence(self, mock_parse_args, mock_exists):
        mock_parse_args.return_value = MockArgs(
            bam1='bam1_file',
            bam2='bam2_file',
            ref1='ref1_file',
            ref2='ref2_file',
            snp='snp_file',
            outprefix='outprefix'
        )

        # 模拟文件不存在
        mock_exists.return_value = False

        with self.assertRaises(SystemExit):
            get_args_and_check_file()

if __name__ == '__main__':
    unittest.main()

