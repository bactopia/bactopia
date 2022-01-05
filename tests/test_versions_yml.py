"""
This script has been adapted from nf-core/modules. Thank you nf-core Team!

Original Source: https://github.com/nf-core/modules/tree/master/tests

"""
from pathlib import Path
import pytest
import yaml
import re
from textwrap import dedent

# 'output' should be the preferred sample name, but it's not always possible
ACCEPTED_IDS = ['output', 'GCF_000292685', 'SRX1390609', 'GCF_000017085']

# some processes have versions in subdirs (e.g. database names)
HAS_SUBDIRS = ['ariba_analysis', 'blast', 'call_variants', 'sequence_type']

def _get_workflow_names():
    """Get all names of all workflows which have a test.yml in the tests directory.
    To do so, recursively finds all test.yml files and parses their content.
    """
    here = Path(__file__).parent.parent.resolve()
    pytest_workflow_files = here.glob("**/test.yml")
    for f in pytest_workflow_files:
        test_config = yaml.load(f.read_text(), Loader=yaml.BaseLoader)
        if test_config:
            for workflow in test_config:
                yield workflow["name"]

def validate_versions_yml(versions_yml):
    """

    """
    assert (
        "END_VERSIONS" not in versions_yml
    ), "END_VERSIONS detected in versions.yml. This is a sign of an ill-formatted HEREDOC"

    # Raises an exception if yaml is not valid
    versions = yaml.safe_load(versions_yml)
    assert (
        len(versions) == 1
    ), "The top-level of versions.yml must contain exactly one entry: the process name as dict key"
    software_versions = next(iter(versions.values()))
    assert len(software_versions), "There must be at least one version emitted."
    for tool, version in software_versions.items():
        assert re.match(
            r"^\d+.*", str(version)
        ), f"Version number for {tool} must start with a number. "


@pytest.mark.workflow(*_get_workflow_names())
def test_ensure_valid_version_yml(workflow_dir):
    workflow_dir = Path(workflow_dir)
    process_name = workflow_dir.name.split("-")[0]
    tested = False
    try:
        if Path(f"{workflow_dir}/bactopia-tools/").exists():
            # Treat as bactopia-tool
            versions_yml_files = Path(workflow_dir / f"bactopia-tools/{process_name}/{process_name}").rglob("*versions.yml")
            for versions_yml_file in versions_yml_files:
                versions_yml = versions_yml_file.read_text()
                validate_versions_yml(versions_yml)
                tested = True
        else:
            # Treat as main workflow
            if process_name in HAS_SUBDIRS:
                for accepted_id in ACCEPTED_IDS:
                    subdirs = Path(workflow_dir / f"pytest/{accepted_id}/logs/{process_name}").glob("*/")
                    for subdir in subdirs:
                        versions_yml_file = subdir / f"versions.yml"
                        versions_yml = versions_yml_file.read_text()
                        validate_versions_yml(versions_yml)
                        tested = True
            else:
                for accepted_id in ACCEPTED_IDS:
                    versions_yml_file = workflow_dir / f"pytest/{accepted_id}/logs/{process_name}/versions.yml"
                    if Path(versions_yml_file).exists():
                        versions_yml = versions_yml_file.read_text()
                        validate_versions_yml(versions_yml)
                        tested = True
                        break
    except FileNotFoundError:
        raise AssertionError(
            dedent(
                f"""\
                `versions.yml` not found in the output directory.
                Expected path: `{versions_yml_file}`

                This can have multiple reasons:
                * The test-workflow failed before a `versions.yml` could be generated.
                * The workflow name in `test.yml` does not start with the tool name.
                """
            )
        )

    assert (
            tested == True
        ), f"No versions.yml was tested. This is a sign the test name is not correctly formatted (got {process_name})"


