import fs from 'fs'
import path from 'path'
import readline from 'readline'
import zlib from 'zlib'

/**
 * Processes a FASTA file and splits it into multiple files based on the part before "::" in the header
 * @param inputFilePath Path to the input FASTA file
 * @param outputDir Directory where output files will be saved
 */
async function processFastaFile(inputFilePath: string, outputDir: string) {
  // Create appropriate stream based on whether file is gzipped
  let fileStream: NodeJS.ReadableStream

  if (inputFilePath.endsWith('.gz')) {
    // For gzipped files, create a gunzip stream
    const gzipStream = fs.createReadStream(inputFilePath)
    const gunzipStream = zlib.createGunzip()
    fileStream = gzipStream.pipe(gunzipStream)
  } else {
    // For regular files, just create a read stream
    fileStream = fs.createReadStream(inputFilePath)
  }

  const rl = readline.createInterface({
    input: fileStream,
    crlfDelay: Infinity, // To handle different line endings
  })

  let currentHeader = ''
  let currentSequence = ''
  let currentOutputFile = ''
  const openFiles = new Set<string>()

  // Process each line in the FASTA file
  for await (const line of rl) {
    if (line.startsWith('>')) {
      // Save the previous sequence if there was one
      if (currentHeader && currentSequence) {
        fs.appendFileSync(
          currentOutputFile,
          `${currentHeader}\n${currentSequence}\n`,
        )
      }

      // Parse the new header
      currentHeader = line
      currentSequence = ''

      // Extract the part before "::"
      const match = line.match(/^>([^:]+)::/)
      if (match && match[1]) {
        const prefix = match[1]
        currentOutputFile = path.join(outputDir, `${prefix}.fa`)

        // Create the file if it doesn't exist yet
        if (!openFiles.has(currentOutputFile)) {
          fs.writeFileSync(currentOutputFile, '', { flag: 'w' })
          openFiles.add(currentOutputFile)
        }
      } else {
        console.warn(`Could not parse header: ${line}`)
      }
    } else if (line.trim() && currentHeader) {
      // Add to the current sequence
      currentSequence += line
    }
  }

  // Save the last sequence
  if (currentHeader && currentSequence && currentOutputFile) {
    fs.appendFileSync(
      currentOutputFile,
      `${currentHeader}\n${currentSequence}\n`,
    )
  }

  console.log(`Processed ${inputFilePath} and created files in ${outputDir}`)
}

const args = process.argv.slice(2)

if (args.length < 1) {
  console.error(
    'Usage: node script.js <input_file_or_directory> [output_directory]',
  )
  process.exit(1)
}

const inputPath = args[0]
const outputDir = args[1] || 'split_fasta_output'

// Create output directory if it doesn't exist
if (!fs.existsSync(outputDir)) {
  fs.mkdirSync(outputDir, { recursive: true })
}

// Process all .fa files in the directory
const files = fs
  .readdirSync(inputPath)
  .filter(file => file.endsWith('.fa.gz'))
  .filter(file => !file.includes('slop'))
  .map(file => path.join(inputPath, file))

for (const file of files) {
  console.log({ file })
  await processFastaFile(file, outputDir)
}
