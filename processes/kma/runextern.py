from logging import Logger
from os import getenv
from pathlib import Path
from subprocess import PIPE, run
from typing import Any, Literal, overload

CONDA_BIN_DIR = getenv("CONDA_BIN_DIR", "/opt/conda/condabin")


def _run_external_text(
    command_line: list[str], logger: Logger, working_dir: Path
) -> tuple[str, str]:
    completed_process = run(
        command_line,
        cwd=working_dir,
        stdout=PIPE,
        stderr=PIPE,
        text=True,
        encoding="utf-8",
    )

#    if completed_process.returncode != 0:
#        error_message = f"Execution of {command_line} failed! Exited with code: {completed_process.returncode}"
#        logger.error(error_message)
#        if completed_process.stderr:
#            for line in completed_process.stderr.split("\n"):
#                logger.error(line)
#        raise RuntimeError(error_message)
    return completed_process.stdout, completed_process.stderr


def _run_external_bytes(
    command_line: list[str], logger: Logger, working_dir: Path
) -> tuple[bytes, bytes]:
    completed_process = run(
        command_line,
        cwd=working_dir,
        stdout=PIPE,
        stderr=PIPE,
        text=False,
    )
    if completed_process.returncode != 0:
        error_message = f"Execution of {command_line} failed! Exited with code: {completed_process.returncode}"
        logger.error(error_message)
        std_err = completed_process.stderr.decode(encoding="utf-8")
        if std_err:
            for line in std_err.split("\n"):
                logger.error(line)
        raise RuntimeError(error_message)
    return completed_process.stdout, completed_process.stderr


@overload
def run_external(
    command_line: list[Any],
    logger: Logger,
    text: Literal[True] = True,
    working_dir: Path = Path(),
) -> tuple[str, str]:
    """
    Run an external command with text outputs
    """


@overload
def run_external(
    command_line: list[Any],
    logger: Logger,
    text: Literal[False],
    working_dir: Path = Path(),
) -> tuple[bytes, bytes]:
    """
    Run an external command with bytes outputs
    """


def run_external(
    command_line: list[Any],
    logger: Logger,
    use_mamba_env: bool = False,
    text: bool = True,
    working_dir: Path = Path(),
) -> tuple[str, str] | tuple[bytes, bytes]:
    """
    Run an external tool
    """
    logger.info("Running external tool %s", command_line)

    if use_mamba_env:
        logger.info(f"Mamba env {CONDA_BIN_DIR} used to run external tool")
        command_line = [f"{CONDA_BIN_DIR}/mamba", "run"] + command_line

    if text:
        std_out, std_err = _run_external_text(
            command_line=command_line, logger=logger, working_dir=working_dir
        )
    else:
        std_out, std_err = _run_external_bytes(
            command_line=command_line, logger=logger, working_dir=working_dir
        )
    logger.info("Done running external tool %s", command_line)
    return std_out, std_err
