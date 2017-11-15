# golgi-polarization

## Plan of attack

1. Import tiff stacks into Python
2. ROI detection/selection
    - Use spray can values with DBSCAN
    - Remove spray can feature   
3. Golgi
    - Issues
        - Large spray can 
            - Solved by discriminatating on radius
        - Beads for purifying/separating endothelial cells
            - Programatically unsolved; will be flagged in a visual inspection.
        - Multiple clusters
            - Half-solved by selecting the cluster closest to wound edge. Sometimes DBSCAN doesn't correctly group Golgi fragments together, and there's a tiny part which is closer to the edge. 
        - Different spray can values at times 2hrs, 4hrs, 6hrs
            - Solved by introducing variable histogram range. 
    - Selection
        - Image segmentation with Otsu's Method
        - DBSCAN clustering
4. Nucleus
    - Issues
        - Overlapping nuclei
            - Solved using selection steps below
    - Selection
        - Binary: Otsu's Method
        - Distance transform
        - Mark maxima
        - Watershed segmentation
5. Export data
    - Convert to a pandas DataFrame
    - Save as csv


### Dependencies
- scikit-image version: 0.13.0
- scipy version: 0.19.1
- matplotlib version: 2.1.0



Microscope: 

Magnification: 

Physical length of a pixel on the CCD

0.108 micron/pixel or 9.23 pixels/micron
