def genomicConditionals() {
    return ["both", "genomic"]
}

def organellarConditionals() {
    return ["both", "organellar"]
}

def isOrganellar(meta) {
    // ORGANELLAR ASSEMBLIES ARE FLAGGED VIA meta.assembly_type
    return meta.assembly_type in ["MITO", "PLASTID"]
}

def runConditionals(meta) {
    // RESOLVE THE run_* CONDITIONAL LIST THAT APPLIES TO THIS ASSEMBLY ITEM
    return isOrganellar(meta) ? organellarConditionals() : genomicConditionals()
}

def getEmptyPlaceholder(idx) {
    // index (idx) as input
    // create an empty placeholder file with index in name
    def placeholder = file("${workDir}/EMPTY_PLACEHOLDER_${idx}")
    if (!placeholder.exists()) {
        placeholder.text = ''
    }
    return placeholder
}
