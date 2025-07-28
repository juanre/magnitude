# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.1.0] - 2025-XX-XX

### Added
- Complete type hints throughout the codebase
- New exception hierarchy for better error handling:
  - `UnitError` - Raised when an unknown unit string is encountered
  - `IncompatibleUnitsError` - Raised for operations on incompatible units, with attributes for accessing the conflicting units and operation
  - `ConversionError` - Raised when a value cannot be converted to a Magnitude
- `__version__` attribute for programmatic version access
- Modern development tooling:
  - Pre-commit hooks with isort, black, ruff, and mypy
  - Pytest configuration to run both unit tests and doctests

### Changed
- **BREAKING**: `mg()` function now raises `IncompatibleUnitsError` when the output unit (`ounit`) has different dimensions than the input unit, instead of silently falling back to the input unit
- **BREAKING**: `Magnitude.coerce()` now raises `ConversionError` instead of returning `None` when it cannot convert a value
- **BREAKING**: Various operations that previously raised generic `MagnitudeError` now raise more specific exception types
- **BREAKING**: Minimum Python version is now 3.9 (was 3.8)
- Migrated to modern Python packaging with `pyproject.toml`
- License changed from Apache 2.0 to MIT (less restrictive)
- Documentation moved from module docstring to README.md
- All string formatting updated from % style to f-strings
- Line length limit increased from 88 to 99 characters (matching Python standard library)

### Removed
- Removed redundant version number from source file header (now only in `pyproject.toml` and `__version__`)

### Fixed
- Fixed `__ipow__` method signature to match `__pow__` for proper type compatibility
- Fixed all type checking issues identified by mypy
- Fixed all reverse arithmetic operations (`__radd__`, `__rsub__`, `__rtruediv__`, etc.) to properly handle non-Magnitude left operands
- Added missing `__rfloordiv__`, `__rmod__`, and fixed `__rdivmod__` implementations
- Added `__repr__` method for better REPL display

## [1.0.2] - Previous version

(Previous changelog entries would go here)
