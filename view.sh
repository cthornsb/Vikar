./dump VIKAR.root VIKAR 3 eject_hitX eject_hitY eject_hitZ
mv dump.out dump2.out
./dump VIKAR.root VIKAR 3 recoil_hitX recoil_hitY recoil_hitZ
vpython ./tools/Layout.py
