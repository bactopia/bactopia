// Import generic module functions
include { initOptions; saveFiles } from '../../../../lib/nf/functions'
options       = initOptions(params.containsKey("options") ? params.options : [:], 'custom_dumpsoftwareversions')
options.btype = "comparative"
options.process_name = "software-versions"
conda_tools   = "bioconda::multiqc=1.25.1"
conda_name    = conda_tools.replace("=", "-").replace(":", "-").replace(" ", "-")
conda_env     = file("${params.condadir}/${conda_name}").exists() ? "${params.condadir}/${conda_name}" : conda_tools

process CUSTOM_DUMPSOFTWAREVERSIONS {
    label 'process_low'

    // Requires `pyyaml` which does not have a dedicated container but is in the MultiQC container
    conda (params.enable_conda ? conda_env : null)
    container "${ workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.25.1--pyhdfd78af_0' :
        'quay.io/biocontainers/multiqc:1.25.1--pyhdfd78af_0' }"

    input:
    path versions

    output:
    path "software_versions.yml", emit: yml
    path "software_versions_mqc.yml", emit: mqc_yml
    path "*.{log,err}", emit: logs, optional: true
    path ".command.*", emit: nf_logs
    path "versions.yml", emit: versions

    script:
    prefix = "software_versions"
    """
    #!/usr/bin/env python
    import datetime
    import yaml
    import platform
    from textwrap import dedent

    def _make_versions_html(versions):
        html = [
            dedent(
                '''\\
                <style>
                #nf-core-versions tbody:nth-child(even) {
                    background-color: #f2f2f2;
                }
                </style>
                <table class="table" style="width:100%" id="nf-core-versions">
                    <thead>
                        <tr>
                            <th> Process Name </th>
                            <th> Software </th>
                            <th> Version  </th>
                        </tr>
                    </thead>
                '''
            )
        ]
        for process, tmp_versions in sorted(versions.items()):
            html.append("<tbody>")
            for i, (tool, version) in enumerate(sorted(tmp_versions.items())):
                html.append(
                    dedent(
                        f'''\\
                        <tr>
                            <td><samp>{process if (i == 0) else ''}</samp></td>
                            <td><samp>{tool}</samp></td>
                            <td><samp>{version}</samp></td>
                        </tr>
                        '''
                    )
                )
            html.append("</tbody>")
        html.append("</table>")
        return "\\n".join(html)

    module_versions = {}
    module_versions["custom_dumpsoftwareversions"] = {
        'python': platform.python_version(),
        'yaml': yaml.__version__
    }

    with open("$versions") as f:
        workflow_versions = yaml.load(f, Loader=yaml.BaseLoader) | module_versions

    workflow_versions["Workflow"] = {
        "Nextflow": "$workflow.nextflow.version",
        "$workflow.manifest.name": "$workflow.manifest.version",
        "command": "$workflow.commandLine",
        "date": datetime.datetime.now()
    }

    versions_mqc = {
        'id': 'software_versions',
        'section_name': '${workflow.manifest.name} Software Versions',
        'section_href': 'https://github.com/${workflow.manifest.name}',
        'plot_type': 'html',
        'description': 'are collected at run time from the software output.',
        'data': _make_versions_html(workflow_versions)
    }

    with open("software_versions.yml", 'w') as f:
        yaml.dump(workflow_versions, f, default_flow_style=False, width=float("inf"))
    with open("software_versions_mqc.yml", 'w') as f:
        yaml.dump(versions_mqc, f, default_flow_style=False, width=float("inf"))

    with open('versions.yml', 'w') as f:
        yaml.dump(module_versions, f, default_flow_style=False, width=float("inf"))
    """
}
