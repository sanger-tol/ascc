def getEmptyPlaceholder(idx) {
    // index (idx) as input
    // create an empty placeholder file with index in name
    def placeholder = file("${workDir}/EMPTY_PLACEHOLDER_${idx}")
    if (!placeholder.exists()) {
        placeholder.text = ''
    }
    return placeholder
}
