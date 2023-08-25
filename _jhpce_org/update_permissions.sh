#!/bin/bash


echo "**** Updating permissions for $1 ****"
date
echo ""

if [[ $HOSTNAME == compute-* ]] || [[ $HOSTNAME == transfer-* ]]; then
    echo "**** Note that warning/error messages are expected for files and directories that you are not the owner of."
    echo "The expected warning/error messages are: "
    echo "    'chgrp: changing group of ‘some_JHPCE_file_path’: Operation not permitted'"
    echo " or 'chmod: changing permissions of ‘some_JHPCE_file_path’: Operation not permitted'."
    echo "If for many files you are not the owner (creator of), you will get lots of these warning/error messages, this is expected!"
    echo "Error or warnings with another syntax are likely real. ****"
    echo ""
    echo "You will need to re-run this script anytime you upload files to JHPCE through Cyberduck / WinSCP as they break the ACLs."
    echo "Every new team member on a given project will likely also need to run this script once."
    echo ""
    echo "For more details about setting permissions at JHPCE using ACLs, please check https://lcolladotor.github.io/bioc_team_ds/organizing-your-work.html#setting-jhpce-file-permissions."
    echo ""
    echo "This message will be displayed for 90 seconds before the script proceeds."
    echo "That way you will have time enough to read it and/or copy it."
    sleep 90
    
    echo ""
    echo "**** Setting read (R), write (W), and execute (X) permissions for hickslab ****"
    sleep 5
    date
    
    find ${1} -user ${USER} -type d -exec nfs4_setfacl -a "A:g:hickslab@cm.cluster:RWX" {} \;
    find ${1} -user ${USER} -type d -exec nfs4_setfacl -a "A:gfdi:hickslab@cm.cluster:RWX" {} \;
    find ${1} -user ${USER} -type f -exec nfs4_setfacl -a "A:g:hickslab@cm.cluster:RW" {} \;
    
    echo ""
    echo "**** Setting read (R), write (W), and execute (X) permissions for lieber_lcolladotor ****"
    sleep 5
    date
    
    find ${1} -user ${USER} -type d -exec nfs4_setfacl -a "A:g:lieber_lcolladotor@cm.cluster:RWX" {} \;
    find ${1} -user ${USER} -type d -exec nfs4_setfacl -a "A:gfdi:lieber_lcolladotor@cm.cluster:RWX" {} \;
    find ${1} -user ${USER} -type f -exec nfs4_setfacl -a "A:g:lieber_lcolladotor@cm.cluster:RW" {} \;
    
    echo ""
    echo "**** Setting read (R), write (W), and execute (X) permissions for lieber_marmaypag ****"
    sleep 5
    date
    
    find ${1} -user ${USER} -type d -exec nfs4_setfacl -a "A:g:lieber_marmaypag@cm.cluster:RWX" {} \;
    find ${1} -user ${USER} -type d -exec nfs4_setfacl -a "A:gfdi:lieber_marmaypag@cm.cluster:RWX" {} \;
    find ${1} -user ${USER} -type f -exec nfs4_setfacl -a "A:g:lieber_marmaypag@cm.cluster:RW" {} \;  
    
    ## To move away from lieber_jaffe
    echo ""
    if getent group lieber_marmaypag | grep -q "\b${USER}\b"; then
        echo "**** Running chgrp lieber_marmaypag ****"
        sleep 5
        date
        chgrp lieber_marmaypag -R ${1}
    elif getent group lieber_lcolladotor | grep -q "\b${USER}\b"; then
        echo "**** Running chgrp lieber_lcolladotor ****"
        sleep 5
        date
        chgrp lieber_lcolladotor -R ${1}
    elif getent group hickslab | grep -q "\b${USER}\b"; then
        echo "**** Running chgrp hickslab ****"
        sleep 5
        date
        chgrp hickslab -R ${1}       
    else
        echo "**** Skipping chgrp step ****"
    fi
    
    ## For setting the group sticky bit
    echo ""
    echo "**** Setting the group sticky bit ****"
    sleep 5
    date
    find ${1} -user ${USER} -type d | xargs chmod g+s
    
    ## Check settings
    echo ""
    echo "**** Checking the nfs4 (ACLs) settings ****"
    sleep 5
    date
    nfs4_getfacl ${1}
else
    echo "**** This script can only work on a qrsh / qsub session. It does not work on a login node since it does not have access to the nfs4_setfacl and nfs4_getfacl commands.****"
fi