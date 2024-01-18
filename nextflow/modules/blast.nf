process BLAST {
    tag "$meta.id"

    input:

    output:

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"    
    """

    """
}

