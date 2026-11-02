/*
================================================================================
    ANNOVAR database maintenance (--mode check_db / update_db, and DB_INSTALL
    before annotation)

    These processes act on the shared humandb directory itself: it is a path
    input (so containers mount it) used through its real path, never a copy.
    DB_DOWNLOAD, DB_CLINVAR and DB_INSTALL write into it, always via a private
    staging directory and an atomic rename after verification.
================================================================================
*/

process DB_AUDIT {
    tag "${build}"
    label 'process_single'

    input:
    tuple val(build), val(entries)
    path humandb
    val stage           // 'before' or 'after' an update; only changes file names

    output:
    tuple val(build), path("${build}.db_audit.${stage}.tsv"), emit: tsv
    path "${build}.db_audit.${stage}.json",                   emit: json
    path "${build}.db_audit.${stage}.txt",                    emit: txt

    script:
    def offline = params.db_offline.toString() == 'true' ? '--offline' : ''
    """
    NO_COLOR=1 annovar_db.py audit \\
        --humandb \$(readlink -f ${humandb}) \\
        --build   ${build} \\
        --entries '${entries}' \\
        --tsv     ${build}.db_audit.${stage}.tsv \\
        --json    ${build}.db_audit.${stage}.json \\
        ${offline} \\
        > ${build}.db_audit.${stage}.txt
    """

    stub:
    """
    printf 'build\\tentry\\tdatabase\\toperation\\tinstalled\\tsize_gb\\tintegrity\\tnewest_official\\tstatus\\taction\\tpackage\\tnotes\\n' > ${build}.db_audit.${stage}.tsv
    echo '{}' > ${build}.db_audit.${stage}.json
    touch ${build}.db_audit.${stage}.txt
    """
}

process DB_DOWNLOAD {
    tag "${build}:${database}"
    label 'process_single'
    maxForks params.db_download_forks.toString().toInteger()

    input:
    tuple val(build), val(database), val(force)
    path humandb
    path annovar_src

    output:
    tuple val(build), val(database), emit: done
    path "${build}_${database}.download.log", emit: log

    script:
    def replace = force ? '--force' : ''
    """
    annovar_db.py download \\
        --humandb     \$(readlink -f ${humandb}) \\
        --build       ${build} \\
        --annovar-src \$(readlink -f ${annovar_src}) \\
        ${replace} \\
        '${database}' 2>&1 | tee ${build}_${database}.download.log
    """

    stub:
    """
    touch ${build}_${database}.download.log
    """
}

process DB_CLINVAR {
    tag "${build}:clinvar_${release}"
    label 'process_medium'

    input:
    tuple val(build), val(release), val(force)    // release: YYYYMMDD or 'latest'
    path humandb

    output:
    tuple val(build), env('CLINVAR_DB'), emit: done
    path "${build}_clinvar.build.log",   emit: log

    script:
    def rebuild = force ? '--force' : ''
    """
    CLINVAR_DB=\$(annovar_db.py clinvar \\
        --humandb \$(readlink -f ${humandb}) \\
        --build   ${build} \\
        --date    ${release} \\
        ${rebuild} \\
        2> >(tee ${build}_clinvar.build.log >&2))
    """

    stub:
    """
    CLINVAR_DB=clinvar_stub
    touch ${build}_clinvar.build.log
    """
}

// Before annotation (--auto_install_db): install the protocol databases of one
// build that are not installed yet; installed ones are left as they are.
/*
 * The report's gene table (humandb/annovar_nf/genes.tsv: GenCC, gnomAD constraint,
 * HGNC) - built when missing, or older than max_age days ('' keeps any table).
 * A failed download is a warning: the report is then made without gene information.
 */
process DB_GENES {
    label 'process_single'

    input:
    path humandb
    val  max_age

    output:
    val true,                        emit: done
    path "genes.db_install.log",     emit: log

    script:
    def age = max_age ? "--max-age ${max_age}" : ''
    """
    annovar_db.py genes --humandb \$(readlink -f ${humandb}) ${age} --optional 2>&1 | tee genes.db_install.log
    """

    stub:
    """
    touch genes.db_install.log
    """
}

process DB_INSTALL {
    tag "${build}"
    label 'process_medium'

    input:
    tuple val(build), val(entries)
    path humandb
    path annovar_src

    output:
    val build,                        emit: build
    path "${build}.db_install.log",   emit: log

    script:
    """
    annovar_db.py install \\
        --humandb         \$(readlink -f ${humandb}) \\
        --build           ${build} \\
        --entries         '${entries}' \\
        --annovar-src     \$(readlink -f ${annovar_src}) \\
        --clinvar-release ${params.clinvar_release} \\
        --forks           ${params.db_download_forks} \\
        2>&1 | tee ${build}.db_install.log
    """

    stub:
    """
    touch ${build}.db_install.log
    """
}
