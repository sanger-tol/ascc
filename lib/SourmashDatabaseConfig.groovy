/**
 * Class for parsing and validating Sourmash database configurations
 *
 * This class is responsible for:
 * 1. Parsing database configurations from CSV files or params.sourmash_databases
 * 2. Validating database configurations (file existence, parameter validity)
 * 3. Logging database configuration summary
 *
 *
 * @see subworkflows/local/run_sourmash/main.nf for runtime logic
 */
class SourmashDatabaseConfig {

    /**
     * Parse database configuration from either CSV file or params.sourmash_databases
     * Priority: params.sourmash_db_config (CSV) > params.sourmash_databases (List of Maps)
     *
     * @param params Pipeline parameters object
     * @return List of Maps with database configurations
     *         Each map contains: [name, path, k_available, k_for_search, s, assembly_taxa_db]
     */
    static List<Map> parseDatabaseConfig(params) {
        // Check if CSV config file is provided
        if (params.sourmash_db_config) {
            return parseDatabasesFromCSV(params.sourmash_db_config)
        }

        // Use params.sourmash_databases if provided
        if (params.sourmash_databases && params.sourmash_databases.size() > 0) {
            return params.sourmash_databases
        }

        // No configuration provided
        return []
    }

    /**
     * Parse database configuration from CSV file
     * Expected CSV format: name,path,k_available,k_for_search,s,assembly_taxa_db
     * k_available should be in format "[21,31,51]"
     *
     * @param csvPath Path to CSV configuration file
     * @return List of Maps with database configurations
     */
    static List<Map> parseDatabasesFromCSV(csvPath) {
        def databases = []
        def csvFile = new File(csvPath)

        if (!csvFile.exists()) {
            System.err.println "[ASCC Sourmash] CSV configuration file not found: ${csvPath}"
            return []
        }

        csvFile.eachLine { line, lineNum ->
            // Skip header line
            if (lineNum == 0) return

            def fields = line.split(',')
            if (fields.size() != 6) {
                System.err.println "[ASCC Sourmash] Skipping malformed line ${lineNum + 1} in CSV: ${line}"
                return
            }

            def name = fields[0].trim()
            def path = fields[1].trim()
            def k_available_str = fields[2].trim()
            def k_for_search = fields[3].trim().toInteger()
            def s = fields[4].trim().toInteger()
            def assembly_taxa_db = fields[5].trim()

            // Parse k_available from "[21,31,51]" format
            def k_available = parseKArrayFromString(k_available_str)

            databases.add([
                name: name,
                path: path,
                k_available: k_available,
                k_for_search: k_for_search,
                s: s,
                assembly_taxa_db: assembly_taxa_db
            ])
        }

        return databases
    }

    /**
     * Parse k array from string format "[21,31,51]"
     *
     * @param kString String representation of k array
     * @return List of integers
     */
    static List<Integer> parseKArrayFromString(String kString) {
        // Remove brackets and split by comma
        def cleanStr = kString.replaceAll(/[\[\]]/, '').trim()
        return cleanStr.split(',').collect { it.trim().toInteger() }
    }

    /**
     * Validate database configurations
     * Checks: file existence, valid parameters, duplicate names
     *
     * @param databases List of database configurations
     * @return Map with validation results: [valid: boolean, warnings: List<String>]
     */
    static Map validateDatabases(List<Map> databases) {
        def warnings = []
        def valid = true

        if (databases.size() == 0) {
            return [valid: false, warnings: ["No Sourmash databases configured"]]
        }

        // Check for duplicate database names
        def names = databases.collect { it.name }
        def duplicates = names.findAll { names.count(it) > 1 }.unique()
        if (duplicates.size() > 0) {
            warnings.add("Duplicate database names found: ${duplicates.join(', ')}")
        }

        // Check each database configuration
        databases.each { db ->
            // Check if database file exists
            def dbFile = new File(db.path)
            if (!dbFile.exists()) {
                warnings.add("Database file not found: ${db.name} at ${db.path}")
                valid = false
            }

            // Check if assembly_taxa_db file exists
            def taxaFile = new File(db.assembly_taxa_db)
            if (!taxaFile.exists()) {
                warnings.add("Assembly taxa DB file not found for ${db.name}: ${db.assembly_taxa_db}")
                valid = false
            }

            // Validate k_for_search is in k_available
            if (!db.k_available.contains(db.k_for_search)) {
                warnings.add("k_for_search (${db.k_for_search}) not in k_available ${db.k_available} for database ${db.name}")
            }

            // Validate scaled parameter
            if (db.s <= 0) {
                warnings.add("Invalid scaled parameter (${db.s}) for database ${db.name}")
            }
        }

        return [valid: valid, warnings: warnings]
    }

    /**
     * Get database configuration summary as formatted string
     * Returns formatted information about loaded databases
     *
     * @param databases List of database configurations
     * @return Formatted string with database information
     */
    static String getDatabaseSummary(List<Map> databases) {
        def summary = "[ASCC Sourmash] Loaded ${databases.size()} database configuration(s):\n"
        databases.each { db ->
            summary += "  - ${db.name}: k=${db.k_for_search}, s=${db.s}\n"
            summary += "    DB: ${db.path}\n"
            summary += "    Taxa: ${db.assembly_taxa_db}\n"
        }
        return summary
    }
}
