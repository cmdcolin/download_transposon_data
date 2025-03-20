import fs from 'fs'
import path from 'path'
import { execSync } from 'child_process'

// Remove and recreate extracts directory
try {
  fs.rmSync('extracts', { recursive: true, force: true })
} catch (err) {
  // Directory might not exist, that's fine
}
fs.mkdirSync('extracts', { recursive: true })

for (const fastaFile of fs
  .readdirSync('DNA_transposon_fasta')
  .filter(file => file.endsWith('.fa.gz'))) {
  console.log(fastaFile)

  const parts2 = fastaFile.split('_')
  const X = path.basename(parts2[1], '.dna')

  // Read the .fai index file
  const faiFile = path.join('DNA_transposon_fasta', `${fastaFile}.fai`)

  try {
    const faiContent = fs.readFileSync(faiFile, 'utf8')
    const lines = faiContent.split('\n').filter(line => line.trim())

    for (const line of lines) {
      const [r] = line.split(/\s+/)

      // Parse parts (equivalent to IFS='::' read -a parts <<<"$r")
      const parts = r.split('::')

      // Execute samtools command
      const cmd = `samtools faidx ${path.join('DNA_transposon_fasta', fastaFile)} ${r} --mark-strand custom," ${parts2[0]}_${X}"`
      const output = execSync(cmd).toString()

      // Append to output file
      fs.appendFileSync(path.join('extracts', `${parts[0]}.fa`), output)
    }
  } catch (err) {
    console.error(`Error processing ${fastaFile}: ${err.message}`)
  }
}
