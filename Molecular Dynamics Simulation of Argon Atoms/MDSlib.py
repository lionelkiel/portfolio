import numpy as np
import matplotlib.pyplot as plt

def potential(distances):
    """Returns numpy array of dimensionless potential energy given array of distances of target atom to all other atoms."""
    potential = np.copy(distances)
    potential[potential!=0] = 4*(potential[potential!=0]**-12 - potential[potential!=0]**-6)
    return(potential)

def force_prefactor(distances):
    """Returns numpy array of dimensionless prefactor of the force given array of distances of target atom to all other atoms."""
    prefactor = np.copy(distances)
    prefactor[prefactor!=0] = 24*(2*prefactor[prefactor!=0]**-14 - prefactor[prefactor!=0]**-8)
    return(prefactor)

def FCClattice(positions, length, a):
    """Sets the initial positions of our atoms as a 3 dimensional FCC lattice.   

    Parameters
    ----------
    postions: The initialised position array
    length: Length of each side the cube in which to fill the FCC lattice configuration
    a: The number of atoms has to be given by 4*a^3 with a an integer. This is because in a cube the number of unit cells of an FCClattice has to be a^3 and they each contain 4 atoms.
    
    """
    #FCC lattice configuration (3 dimensions)
    LatticePoints1 = np.array([(i*length/a,j*length/a,k*length/a) for i in range(a) for j in range(a) for k in range(a)])
    LatticePoints2 = np.array([(i*length/a+length/(2*a),j*length/a+length/(2*a),k*length/a) for i in range(a) for j in range(a) for k in range(a)])
    LatticePoints3 = np.array([(i*length/a+length/(2*a),j*length/a,k*length/a+length/(2*a)) for i in range(a) for j in range(a) for k in range(a)])
    LatticePoints4 = np.array([(i*length/a,j*length/a+length/(2*a),k*length/a+length/(2*a)) for i in range(a) for j in range(a) for k in range(a)])
    LatticePoints = np.concatenate((LatticePoints1,LatticePoints2,LatticePoints3,LatticePoints4))
    
    #Setting the initial positions
    positions[:,:,0] = LatticePoints

def run_simulation(time_steps, Position, Velocity, Force, time_step_length, length):
    """Run the simulation for the given time_steps

    Parameters
    ----------
    time_steps: total amount of time steps to run the simulation for
    Position: numpy array to fill with the position of each atom at each time step
    Velocity: numpy array to fill with the velocity of each atom at each time step
    Force: numpy array to fill with the total force on each atom at each time step
    time_step_length: Length of each time step in unit time
    length: Length of rach side of the box
    
    """
    for t in range(time_steps-1):
            Position[:,:,t+1] = (Position[:,:,t]+time_step_length*Velocity[:,:,t]+time_step_length**2*Force[:,:,t]/2)%length #Calculate position of next step, modulo L so the particle stays in the window
            distance_vectors = (Position[:,:,t+1]- Position[:,:,t+1][:,np.newaxis]+length/2)%length - length/2 #array of distance vectors between each atom taking into account the nearest window principle
            prefactor = force_prefactor(np.linalg.norm(distance_vectors,axis=2))
            force_contributions = prefactor[:,:,np.newaxis]*distance_vectors #Gives an array of the force contributions of each atom for the respective target atom
            Force[:,:,t+1] = np.sum(force_contributions, axis=0)
            Velocity[:,:,t+1] = Velocity[:,:,t] + time_step_length*(Force[:,:,t+1]+Force[:,:,t])/2 #Calculate velocity of next time step
            
    
def Lambda_Algorithm(test_time_steps, Position, Velocity, Force, time_step_length, length, Theoretical_Kinetic_Energy, precision = 0.01):
    """Runs the simulation for a shorter time (given by test_time_steps) to get initial conditions in equilibrium. See documentaion for run_simulation.
    Lambda algorithm will stop when Lambda is 1 +/- precision"""
    Lambda = 10
    while (np.abs(Lambda-1)>precision):
        run_simulation(test_time_steps,Position, Velocity, Force, time_step_length, length)
        
        #Calculate total dimensionless kinetic energy at last time step
        Experimental_Kinetic_Energy = 1/2*np.sum(np.linalg.norm(Velocity[:,:,test_time_steps-1],axis=1)**2)
        Lambda = np.sqrt(Theoretical_Kinetic_Energy/Experimental_Kinetic_Energy)
        print('Value of Lambda: ', Lambda)
        
        #Initializing the new posistions, velocities and forces
        Position[:,:,0] = Position[:,:,test_time_steps-1]
        Velocity[:,:,0] = Velocity[:,:,test_time_steps-1]*Lambda
        Force[:,:,0] = Force[:,:,test_time_steps-1]

def simulate(number, simulation_time, time_step_length, temperature, density, dimension = 3):
    """Simulate the positions and velocities of atoms in a box with a periodic boundary.
    
    Parameters
    ----------
    number: The number of atoms in our simulation, should be of the form 4*a^3 so that it can be initialised as an FCClattice
    simulation_time: Simulation time in dimensionless time
    timestep_length: Timestep length in dimensionless time
    temperature: Dimensionless temperature of simulation
    density: Dimensionless density of simulation
    dimension: Currently we can only initialize our FCClattice in 3 dimensions therefore 3 is set as the default
    
    Returns
    -------
    Position: numpy array of position of each atom at each time step
    Velocity: numpy array of velocitie of each atom at each time step
    Force: numpy array of the total force on each atom at each time step
    Kinetic_Energy: numpy array of the total kinetic energy in each time step
    Potential_Energy: numpy array of the total potential energy in each time step
    Pair_correlation: array containing number of atoms in each of the bins
    Bins: array of bins for the pair correlation function
    Distances: numpy array of each distance at each time step
    Pressure: numpy array of the pressure at each time step
    
    """
    print("Initialization...")
    
    time_steps = int(simulation_time/time_step_length) #Number of timesteps
    length = (number/density)**(1/dimension) #Dimensionless length of our box
    
    #Initializing positions and velocities, forces and energies
    Position = np.full((number,dimension,time_steps),np.nan)
    Velocity = np.full((number,dimension,time_steps),np.nan)
    Force = np.full((number,dimension,time_steps), np.nan)
    
    #Calculate a and return an error if a is not of the form 4*a^3
    a = np.cbrt(number/4)
    if a.is_integer()==False:
        print('ERROR: number needs to be of the form 4*a^3')
        return
    a = int(a)
    
    FCClattice(Position, length, a) #Set initial positions as FCClattice
    Velocity[:,:,0] = np.random.normal(0,np.sqrt(temperature),(number,dimension)) #Maxwell distributed initial velocities
    
    #Calculating initial forces by broadcasting
    distance_vectors = (Position[:,:,0]- Position[:,:,0][:,np.newaxis]+length/2)%length - length/2 #array of distance vectors between each atom taking into account the nearest window principle
    prefactor = force_prefactor(np.linalg.norm(distance_vectors,axis=2))
    force_contributions = prefactor[:,:,np.newaxis]*distance_vectors #Gives an array of the force contributions of each atom for the respective target atom
    Force[:,:,0] = np.sum(force_contributions,axis=0)
    
    print("Performing Lambda algorithm...")
    Theoretical_Kinetic_Energy = (number-1)*3*temperature/2
    test_time_steps = int(time_steps/5) #we take a fifth of our time steps for the lambda algorithm
    Lambda_Algorithm(test_time_steps, Position, Velocity, Force, time_step_length, length, Theoretical_Kinetic_Energy)
    
    print("Performing Simulation...")
    run_simulation(time_steps, Position, Velocity, Force, time_step_length, length)
    
    #Calculate kinetic and potential energy.
    Kinetic_Energy = 1/2*np.sum(np.linalg.norm(Velocity,axis=1)**2,axis=0)
    flattened_Position = Position.reshape(Position.shape[0],Position.shape[1]*Position.shape[2]) #Put time steps and x,y,z axis into one so that the subtraction does what we want
    distance_vectors = ((flattened_Position - flattened_Position[:,np.newaxis]+length/2)%length - length/2).reshape(Position.shape[0],Position.shape[0],Position.shape[1],Position.shape[2]) #Split time steps and x,y,z again
    Distances = np.linalg.norm(distance_vectors, axis = 2)
    Potential_Energy = np.sum(potential(Distances),axis = (0,1))/2
    
    #Calculating Histogram
    Max_bin_edge = np.sqrt(3)*length/2 #The maximum length with the nearest windows principle
    Number_of_bins = 1000
    Length_of_bin = Max_bin_edge/Number_of_bins
    Bins = [i*Length_of_bin for i in range(Number_of_bins+1)] #+1 For edge of last bin
    Histogram,Bins = np.histogram(Distances[Distances!=0], Bins) #We do not wan't to include the distances of atoms to themselves
    Histogram = Histogram/(2*time_steps) #Averaging over the time steps. Also each distance is calculated twice
    Pair_correlation = (Histogram*2*length**3)/(number*(number-1)*4*np.pi*Length_of_bin*Bins[1:]**2) #Take r as right edge of bins, for lots of bins it doesn't matter wheter r is the center or the edge of a bin. For convenience and not having to deal with 0 we take the right edge
    
    #Calculating Pressure
    Pressure = density*(temperature + np.sum(force_prefactor(Distances)*Distances**2,axis = (0,1))/(12*number))
    
    print("Simulation Complete!")
    
    return Position,Velocity,Force,Kinetic_Energy,Potential_Energy,Pair_correlation,Bins,Distances,Pressure
    