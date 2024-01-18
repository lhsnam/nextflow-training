process MAKE_BLAST_DATABASE {

    input:
    file database

    output:
    each '$database*'
    path('*.log')
    path "versions.yml", emit: versions

    script: 
    """
    makeblastdb -in $database -dbtype nucl > make_database.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        makblastdb: \$(makeblastdb -version | sed 's/makeblastdb //g')
    END_VERSIONS
    """
}
