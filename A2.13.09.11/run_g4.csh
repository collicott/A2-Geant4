#!/bin/tcsh

# EXPERIMENT DIRECTORY
set experiment = compton
#set experiment = pi0

# PROCESS
#set processes = "pi0"
set processes = "compton pi0"
#set processes = eff

# TARGET
#set targets = "p c"
#set targets = "p g"
#set targets = g
set targets = p
#set targets = c

# TARGET MATERIAL
set tmat = LH2
#set tmat = PolTarg

# TARGET LENGTH
set len = 10cm
#set len = 5cm
#set len = 2cm

if ( $tmat == "PolTarg") then
	set enclosure = $tmat
else
	set enclosure = "$tmat"_"$len"
endif

# ENERGY
set energies = 250
#set energies = "200 210 220 230"

rm -f run_g4.log

rm -f macros/DetectorSetup.mac
ln -s DetSetup_"$enclosure".mac macros/DetectorSetup.mac

foreach proc ($processes)

	foreach tgt ($targets)

		foreach eg ($energies)

			echo Analysing $proc $tgt $eg
			echo Analysing $proc $tgt $eg >> run_g4.log

			echo /A2/physics/Physics QGSP_BIC > batch.mac
			echo /run/initialize >> batch.mac
			echo /A2/generator/Seed 1111111 >> batch.mac

			if ( ( $proc == "compton") || ( $proc == "eff") ) then
				set np = 1
			else
				set np = 2
			endif
			if ( ( $tgt == "p") && ( $proc != "eff") ) then
				set np = `expr $np + 1`
			endif

			echo /A2/generator/NToBeTracked $np >> batch.mac
			if ( ( $proc == "compton") || ( $proc == "eff") )  then
				echo /A2/generator/Track 1 >> batch.mac
				if  ( ( $tgt == "p") && ( $proc != "eff") ) then
					echo /A2/generator/Track 2 >> batch.mac
				endif
			else
				if  ( $tgt == "p") then
					echo /A2/generator/Track 1 >> batch.mac
				endif
				echo /A2/generator/Track 3 >> batch.mac
				echo /A2/generator/Track 4 >> batch.mac
			endif

			# input file
			set indir = evgen/$len
			echo /A2/generator/InputFile $indir/"$proc"_"$tgt"_"$eg"_in.root \
				>> batch.mac

			# output file
			set outdir = out/$experiment/$enclosure
			echo /A2/event/setOutputFile $outdir/"$proc"_"$tgt"_"$eg"_out.root \
				>> batch.mac

			time A2 batch.mac >>&! run_g4.log &

			wait
		end
	end
end

#time A2 nicegen_comp.mac >>&! run_g4.log &
#wait
#time A2 nicegen_pi0.mac >>&! run_g4.log &
#wait

exit 0
