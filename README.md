# pipeline4CaImg2dFF0
A Pipeline for Calcium Imaging Data to dFF0 Data (via NoRMCorre and suite2p)

The version of MATLAB should be higher than 2022a. (`LinearExtrap.m` used `@pagemldivide`)


## 2023/05/15 ---update---
* change `...\Functions\NoRMCorre-DIY\Additional_Functions\read_bin_file_DIY.m` and `motCRCTviaNoRMCorre.m`
* add `dFF0extract.m` and change `...\Functions\Preprocessing\reExtractF.m`: could re-extract F from data.bin
* change `checkROIdtctdVIAsuite2p.m`: repair function `XY2bwI`
* change `...\suite2p_files_debugged_and_DIY_v0.12.1[@ZetaEta]\gui\merge.py`: merged lambda were corrected
* add `...\suite2p_files_debugged_and_DIY_v0.12.1[@ZetaEta]\logo` (could use `logo.ico` for your suite2p Launcher)
