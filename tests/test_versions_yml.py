"""
This script has been adapted from nf-core/modules. Thank you nf-core Team!

Original Source: https://github.com/nf-core/modules/tree/master/tests

"""
from pathlib import Path
import pytest
import yaml
import re
from textwrap import dedent

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

@pytest.mark.workflow(*_get_workflow_names())
def test_ensure_valid_version_yml(workflow_dir):
    workflow_dir = Path(workflow_dir)
    process_name = workflow_dir.name.split("-")[0]
    try:
        versions_yml_file = workflow_dir / f"output/logs/{process_name}/versions.yml"
        versions_yml = versions_yml_file.read_text()
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
