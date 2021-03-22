$argv1Fichero=@("Usuario-1.txt","Usuario-2.txt","Usuario-3.txt","Usuario-4.txt","Usuario-5.txt","Usuario-6.txt")
$argv5Selector=@(2)
$argv2PopSize=@(100)
$argv3PCross=@(0.85)
$argv4PMut=@(0.075)


for($i=0;$i -lt $argv1Fichero.Length; $i++){
	$total=1

	for($j=0;$j -lt $argv5Selector.Length; $j++){

		for($k=0;$k -lt $argv2PopSize.Length;$k++){

			for($l=0;$l -lt $argv3PCross.Length;$l++){

				for($x=0;$x -lt $argv4PMut.Length;$x++){

					"Ejecutando caso " + $total + " de 80 para fichero " + $argv1Fichero[$i]
					& "C:\Users\diego\Desktop\P1\bin\Debug\P1" $argv1Fichero[$i] $argv2PopSize[$k] $argv3PCross[$l] $argv4PMut[$x] $argv5Selector[$j]
					$total++
				}
			}
		}
	}
}