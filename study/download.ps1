$sourcePath = "~/rangeexp/resk/data/ooa_move3.re"
$destFolder = $PSScriptRoot + "\data\"
$compName = "***.riken.jp"
$cred = "***"
$privateKeyPath = "C:\Users\***\Documents\***\private"

Get-SCPItem -ComputerName $compName -Credential $cred -KeyFile $privateKeyPath -PathType File -Destination $destFolder -Path $sourcePath
