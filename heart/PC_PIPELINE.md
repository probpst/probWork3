# Heart case pipeline on the PC (frames -> probCFD case)

    W=$HOME/heartCFD/work; mkdir -p $W/hg $W/cyc $W/frames
    # scripts with their hard-coded /home/claude paths redirected to $W
    for f in heart/scripts/grid/*.py heart/scripts/restore_state.py; do sed "s#/home/claude#$W#g" $f > $W/hg/$(basename $f); done
    for f in heart/scripts/cycle/*.py;  do sed "s#/home/claude#$W#g" $f > $W/cyc/$(basename $f); done
    pip3 install --user pyamg            # (numpy, scipy already there)
    python3 $W/hg/restore_state.py heart/data          # heart_grid.g + .bcs + block_list -> $W/hg/assembled_s5.pkl (~1 min)
    cd $W/cyc && GRID=$W/hg/assembled_s5.pkl OUT=$W/cyc/bmesh.npz MAP=$W/cyc/vertex_to_grid_nodes.txt python3 bmesh.py
    cmp $W/cyc/vertex_to_grid_nodes.txt heart/data/vertex_to_grid_nodes.txt && echo "vertex order identical"
    python3 tools/morph_wall_frames.py --hg $W/hg --cyc $W/cyc --bmesh $W/cyc/bmesh.npz --wall heart/data/wall_frames --out $W/frames --workers 6
    # (ends with the sliver fix of tools/sliver_fix.py: nodes joined by reference edges < 0.1 mm get one displacement;
    #  raw frames stay in $W/frames/raw; expect 'ALL FRAMES VALID', min scaled Jacobian about 0.036)
    python3 tools/make_heart_case.py cases/heart --grid heart/data/heart_grid.g --bcs heart/data/heart_grid.bcs \
            --frames $W/frames --numbering $W/hg/node_numbering.npz --csv heart/data/cycle_frames.csv --nranks 20
    cp heart/data/valves/*.dat cases/heart/        # valve zones (the $Valves group is added to the namelist separately)
