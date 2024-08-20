$sourcePath = "~/rangeexp/resk/data/ooa_move3.re"
$destFolder = $PSScriptRoot + "\data\"
$compName = "raiden.riken.jp"
$cred = "hartree"
$privateKeyPath = "C:\Users\Hartree\Documents\raiden\private"

Get-SCPItem -ComputerName $compName -Credential $cred -KeyFile $privateKeyPath -PathType File -Destination $destFolder -Path $sourcePath
