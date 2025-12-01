def getEmptyPlaceholder(idx) {
    def placeholder = file("${workDir}/EMPTY_PLACEHOLDER_${idx}")
    if (!placeholder.exists()) {
        placeholder.text = ''
    }
    return placeholder
}
