folderlist=`ls *.xml.gz`

for folder in $folderlist;
  do 
    folder=${folder/.xml.gz/}
    echo $folder
    if [[ $folder == *"HLV"* ]]; then 
      echo ifos = [\'H1\', \'L1\', \'V1\'];
    else echo ifos = [\'H1\', \'L1\']
    fi
    echo webdir = /home/abhirup.ghosh/public_html/testGR_IR/runs/systematics_error_characterisation/${folder}
    echo baseurl = https://ldas-jobs.ligo.caltech.edu/~abhirup.ghosh/testGR_IR/runs/systematics_error_characterisation/${folder}
    echo channels = \{\'H1\':\'H1:${folder}_H-H1HWINJ\',\'L1\':\'L1:${folder}_L-L1HWINJ\',\'V1\':\'V1:${folder}_V-V1HWINJ\'\}
    echo H1-psd = /home/abhirup.ghosh/Documents/Work/testGR_IR/runs/systematics_error_characterisation/${folder}/${folder}_H1_psd.txt
    echo L1-psd = /home/abhirup.ghosh/Documents/Work/testGR_IR/runs/systematics_error_characterisation/${folder}/${folder}_L1_psd.txt
    echo V1-psd = /home/abhirup.ghosh/Documents/Work/testGR_IR/runs/systematics_error_characterisation/${folder}/${folder}_V1_psd.txt
  done

