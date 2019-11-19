#!/usr/bin/env python3

import pytest
from serotools import cli


def test_error_on_empty_command_line():
    """Verify exception on empty command line."""
    with pytest.raises(SystemExit):
        cli.run_from_line("")
