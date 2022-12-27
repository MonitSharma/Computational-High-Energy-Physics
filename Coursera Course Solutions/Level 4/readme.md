# Week4 assignment
In this task you have to design a prediction model that is able to discriminate the basetracks beloning to the electromagnetic showers from the background basetracks. The goal is to get the best possible score (ROC AUC).

## The data files:
At Week4 Release of this repository you can find the following files:

https://github.com/hse-aml/hadron-collider-machine-learning/releases/download/Week_4/training.tgz - archive with brick volumes filled with labeled basetracks https://github.com/hse-aml/hadron-collider-machine-learning/releases/download/Week_4/test.h5.gz - volume with basetracks you have to discriminate

1. Follow the details at week 3 notebook

2. Each BaseTrack (BT) is described by:

3. Coordinates: (X, Y, Z)

4. Angles in the brick-frame: (TX, TY)

5. After the 'add_neighbours' we add dX,dY,dZ,dTX,dTY (distance from the origin) and group close tracks from neighbour plates into pairs.
