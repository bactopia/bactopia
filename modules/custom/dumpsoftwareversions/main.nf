process CUSTOM_DUMPSOFTWAREVERSIONS {
    label 'process_low'

    // Requires `pyyaml` which does not have a dedicated container but is in the MultiQC container
    conda "${task.ext.env.condaDir}/${task.ext.env.toolName}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.env.image : task.ext.env.docker }"

    input:
    path versions

    output:
    path "software_versions.yml"           , emit: yml
    path "software_versions_mqc.yml"       , emit: mqc_yml
    tuple val(meta), path("*.{log,err}")   , emit: logs, optional: true
    tuple val(meta), path(".command.begin"), emit: nf_begin
    tuple val(meta), path(".command.err")  , emit: nf_err
    tuple val(meta), path(".command.log")  , emit: nf_log
    tuple val(meta), path(".command.out")  , emit: nf_out
    tuple val(meta), path(".command.run")  , emit: nf_run
    tuple val(meta), path(".command.sh")   , emit: nf_sh
    tuple val(meta), path(".command.trace"), emit: nf_trace
    tuple val(meta), path("versions.yml")  , emit: versions

    script:
    prefix = "software_versions"

    // Create a new meta variable
    meta = [:]
    meta.id = "${prefix}-${task.process}"
    meta.name = prefix
    meta.output_dir = "${prefix}/tools/${task.ext.process_name}/${task.ext.subdir}"
    meta.logs_dir = "${prefix}/tools/${task.ext.process_name}/${task.ext.subdir}/logs/${task.ext.logs_subdir}"
    meta.process_name = task.ext.process_name
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
