# write function that calculates the acceleration of particle 1 in the x direction, given the input of the the coordinates of two particles: r1 and r2, which are each an array of two items.
    # represent x1 as r1[0] and x2 as r2[0]; represent y1 as r1[1] and x2 as r2[1].
    # return the result of the calculation: -24(x2 - x1)/((x2 - x1)^2 + (y2 - y1)^2)^7 + 12(x2 - x1)/((x2 - x1)^2 + (y2 - y1)^2)^4
# write another function that outputs the acceleration of particle 2 in the x direction, which is the negative of the previously written function.
# write two more functions that get the acceleration for particle 1 and 2 in the y direction, by replacing x2 - x1 in the numerators (of the x-dimensional acceleration calculation) to y2 - y1.

# create an empty array for time.
# create 4 arrays for velocities:
    # an array for velocity of particle 1 in the x direction, with 0 as the first and only term for now,
    # an array for velocity of particle 2 in the x direction, with 0 as the first and only term for now,
    # an array for velocity of particle 1 in the y direction, with 0 as the first and only term for now,
    # an array for velocity of particle 2 in the y direction, with 0 as the first and only term for now.
# create two arrays which are the initial coordinates of r1 and r2.
# create 4 arrays for positions:
    # an array for the x position of particle 1, with r1[0] as the first and only term for now,
    # an array for the y position of particle 1, with r1[1] as the first and only term for now,
    # an array for the x position of particle 2, with r1[0] as the first and only term for now,
    # an array for the y position of particle 2, with r1[1] as the first and only term for now.

# iterate through the range of the number of steps, which is 100.
    # if it is the first iteration, hence index is 0,
        # define 4 new variables representing the updated velocities:
            # particle 1's x-dimensional velocity is the initial velocity 0 plus 1/2 * 0.01 * (particle 1's acceleration in the x direction, given previously defined initial positions of particle 1 and 2),
            # particle 1's y-dimensional velocity is the initial velocity 0 plus 1/2 * 0.01 * (particle 1's acceleration in the y direction, given previously defined initial positions of particle 1 and 2),
            # particle 2's x-dimensional velocity is the initial velocity 0 plus 1/2 * 0.01 * (particle 2's acceleration in the x direction, given previously defined initial positions of particle 1 and 2),
            # particle 2's y-dimensional velocity is the initial velocity 0 plus 1/2 * 0.01 * (particle 2's acceleration in the y direction, given previously defined initial positions of particle 1 and 2).
        # add the 4 new variables into their respective arrays that were previously defined.
        # append the index of the iteration to the time array.
    # if the iteration index is not 0,
        # update the positions for particle 1 and particle 2 by adding (0.01 * the last term of each velocity array) to their respective index in the r1 and r2 arrays.
        # append the x and y coordinates for particle 1 and particle 2 to the four position arrays, which are previously defined.
        # define 4 new variables representing the k for each particle in each dimension:
            # k1x is (0.01 * particle 1's acceleration in the x direction, given the newly updated positions of particle 1 and 2),
            # k1y is (0.01 * particle 1's acceleration in the y direction, given the newly updated positions of particle 1 and 2),
            # k2x is (0.01 * particle 2's acceleration in the x direction, given the newly updated positions of particle 1 and 2),
            # k2y is (0.01 * particle 2's acceleration in the y direction, given the newly updated positions of particle 1 and 2),
        # define 4 new variables representing the updated velocities:
            # particle 1's x-dimensional velocity is the latest velocity (the last term in the velocity array) plus k1x,
            # particle 1's y-dimensional velocity is the latest velocity (the last term in the velocity array) plus k1y,
            # particle 2's x-dimensional velocity is the latest velocity (the last term in the velocity array) plus k2x,
            # particle 2's y-dimensional velocity is the latest velocity (the last term in the velocity array) plus k2y.
        # add the 4 new velocities into their respective arrays.
        # append the index of the iteration to the time array.

# plot the trajectory with the 4 position arrays, which are now filled.
