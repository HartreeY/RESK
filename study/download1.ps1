$sourcePath = "~/rangeexp/resk/data/2d_octfst.re"
$destFolder = $PSScriptRoot + "\data\"
$compName = "raiden.riken.jp"
$cred = "hartree"
$privateKeyPath = "C:\Users\Hartree\Documents\raiden\private"

Get-SCPItem -ComputerName $compName -Credential $cred -KeyFile $privateKeyPath -PathType File -Destination $destFolder -Path $sourcePath
