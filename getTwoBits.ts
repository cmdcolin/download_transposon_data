import fs from 'node:fs'

const files = JSON.parse(fs.readFileSync('r.json', 'utf8'))

console.log(
  files.list.map(r => r.repeats[0].adapter.bigBedLocation.uri).join('\n'),
)
