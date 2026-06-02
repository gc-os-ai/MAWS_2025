"""
Unit tests for the maws2023 command-line argument parser.

These tests exercise argparse only (no LEaP/OpenMM build).
"""

from __future__ import annotations

import sys

import pytest

from maws.maws2023 import parse_args


def test_cli_default_salt_conc(monkeypatch) -> None:
    """--salt-conc defaults to physiological 0.15 mol/L."""
    monkeypatch.setattr(sys, "argv", ["maws"])
    args = parse_args()
    assert args.salt_conc == 0.15


def test_cli_custom_salt_conc(monkeypatch) -> None:
    """--salt-conc overrides the default."""
    monkeypatch.setattr(sys, "argv", ["maws", "--salt-conc", "0.3"])
    args = parse_args()
    assert args.salt_conc == 0.3


def test_cli_zero_salt_conc(monkeypatch) -> None:
    """--salt-conc 0 (documented unscreened mode) is accepted at the boundary."""
    monkeypatch.setattr(sys, "argv", ["maws", "--salt-conc", "0"])
    args = parse_args()
    assert args.salt_conc == 0.0


def test_cli_rejects_negative_salt_conc(monkeypatch) -> None:
    """Negative --salt-conc is rejected by argparse (non-negative-float validator)."""
    monkeypatch.setattr(sys, "argv", ["maws", "--salt-conc", "-0.1"])
    with pytest.raises(SystemExit):
        parse_args()
