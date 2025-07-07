# FoldX Installation

This directory should contain the FoldX binary.

## Download FoldX

1. Visit https://foldxsuite.crg.eu/
2. Register for an academic license (free for academic use)
3. Download the appropriate version for your system
4. Place the binary in this directory
5. Make it executable: `chmod +x foldx_51binaryMac`

## Expected Structure

After downloading, this directory should contain:
```
foldx/
├── README.md (this file)
├── foldx_51binaryMac      # For macOS
└── molecules/             # FoldX molecules directory
```

For Linux systems, the binary might be named `foldx` or `foldx_5`.

## Usage

The GearBind scripts will automatically look for FoldX in this directory.
The default path is hardcoded to `./foldx/foldx_51binaryMac`.

If your FoldX binary has a different name or location, you can override it:
```bash
python script/generate_mutants_foldx.py \
    --foldx-path /path/to/your/foldx \
    ...other arguments...
``` 