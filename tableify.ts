import fs from "node:fs";
import path from "node:path";

const r = JSON.parse(fs.readFileSync("r.json", "utf8"));
for (const item of r.list) {
  console.log(
    [
      path
        .basename(item.repeats[0].adapter.bigBedLocation.uri)
        .replace(".bb", ".bed"),
      item.assembly.sequence.adapter.uri,
      item.assembly.name,
    ].join("\t"),
  );
}
