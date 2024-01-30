import numpy as np
import matplotlib.pyplot as plt
import os
import time
import multiprocessing as mp
import concurrent.futures as cf
import joblib as job


def calculate_stress_and_gyration_tensors(polymer):
    #Calculate the displacement vector between neighboring beads
    dr = np.diff(polymer, axis=0)
    #Apply periodic boundary conditions
    dr = np.where(dr > 35, dr - 70, dr)
    dr = np.where(dr < -35, dr + 70, dr)
    #Calculate the shear stress tensor
    tau_xy = np.sum(dr[:, 0] * dr[:, 1])
    #Calculate the normal stress tensor
    sigma_xx = np.sum(dr[:, 0] ** 2) / polymer.shape[0]
    sigma_yy = np.sum(dr[:, 1] ** 2) / polymer.shape[0]
    #Calculate the pressure tensor
    pressure_tensor = np.array([[sigma_xx, 0], [0, sigma_yy]])
    #Combine the normal and shear stress tensors into a stress tensor
    stress_tensor = np.array([[sigma_xx, tau_xy], [tau_xy, sigma_yy]])
    #Calculate the gyration tensor
    center_of_mass = np.mean(polymer, axis=0)
    gyration_tensor = np.zeros((2, 2))
    for bead in polymer:
        delta = bead - center_of_mass
        gyration_tensor += np.outer(delta, delta)
    gyration_tensor /= polymer.shape[0]
    
    #Calculate eigenvalues and eigenvectors of the stress tensor
    stress_eigenvals, stress_eigenvecs = np.linalg.eig(stress_tensor)
    
    #Calculate eigenvalues and eigenvectors of the pressure tensor
    pressure_eigenvals, pressure_eigenvecs = np.linalg.eig(pressure_tensor)
    
    #Calculate eigenvalues and eigenvectors of the gyration tensor
    gyration_eigenvals, gyration_eigenvecs = np.linalg.eig(gyration_tensor)
    
    # Calculate the squared radius of gyration
    rg_sq = np.sum(gyration_eigenvals**2) #/ len(polymer)

    # Calculate the radius of gyration
    rg = np.sqrt(rg_sq)
    rg2 = rg.reshape((1, 1)) # Rg reshaped
                    
    # Return the stress tensor, pressure tensor, gyration tensor, eigenvalues, and eigenvectors
    return stress_tensor, pressure_tensor, gyration_tensor, stress_eigenvals, stress_eigenvecs, pressure_eigenvals, pressure_eigenvecs, gyration_eigenvals, gyration_eigenvecs,  rg2 
    #        0                 1                 2              3                    4                5                   6                      7                   8           9

def compute_angles(eigenvec1, eigenvec2):

    #Normalize the eigenvectors
    eigenvec1 /= np.linalg.norm(eigenvec1)
    eigenvec2 /= np.linalg.norm(eigenvec2)

    #Compute the angle between the two eigenvectors
    angle = np.arccos(np.clip(np.dot(eigenvec1, eigenvec2), -1.0, 1.0))

    #Convert angle from radians to degrees
    angle_deg = np.degrees(angle)
    return angle_deg

def get_lock(lock_file, max_wait=60):
    """Acquire a lock by creating a lock file, or wait until the lock is released."""
    lock_file = lock_file.strip('.lock') + '.lock'
    start_time = time.time()
    while os.path.exists(lock_file):
        if time.time() - start_time > max_wait:
            raise TimeoutError(f"Failed to acquire lock on {lock_file}")
        time.sleep(0.1)
    open(lock_file, 'w').close()  #create lock file
    return Lock(lock_file)
    
class Lock:
    """A lock that releases the lock file on exit."""
    def __init__(self, lock_file):
        self.lock_file = lock_file

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        os.remove(self.lock_file)
            
def polymer_analysis_wrapper(args):
    time_frame, coords_all = args
    coords = coords_all[time_frame]
    for polymer in coords:
        results = calculate_stress_and_gyration_tensors(polymer)

        #Compute and write out the angles
        max_eigenvec_stress_idx = np.argmax(np.abs(results[3]))
        max_eigenvec_stress = results[4][:, max_eigenvec_stress_idx]

        max_eigenvec_gyration_idx = np.argmax(np.abs(results[7]))
        max_eigenvec_gyration = results[8][:, max_eigenvec_gyration_idx]

        angle = compute_angles(max_eigenvec_stress, max_eigenvec_gyration)

        #max_eigenvec_stress = results[4][:, np.argmax(results[3])]
        #max_eigenvec_gyration = results[8][:, np.argmax(results[7])]

        #angle = compute_angles(max_eigenvec_stress, max_eigenvec_gyration)

        #compute the Max Eigen values gyration & stress Tensors
        max_index_gyration = np.argmax(np.abs(results[7]))
        max_eigenvalue_gyration = results[7][max_index_gyration]
        sign_gyration = np.sign(results[7][max_index_gyration])
        max_eigenvalue_gyration *= sign_gyration

        max_index_stress = np.argmax(np.abs(results[3]))
        max_eigenvalue_stress = results[3][max_index_stress]
        sign_stress = np.sign(results[3][max_index_stress])
        max_eigenvalue_stress *= sign_stress

        #max_eigenvalue_gyration = np.max(results[7])
        #max_eigenvalue_stress = np.max(results[3])

        #with get_lock(f'Tensor_Analysis/Stress_Tensor_Output.dat.lock'):
        with open(f'Tensor_Analysis/Stress_Tensor_Output.dat', 'a') as fp1:
            np.savetxt(fp1, results[0], delimiter='        ', fmt='%.14F')
 
        #with get_lock(f'Tensor_Analysis/Pressure_Tensor_Output.dat.lock'):
        with open(f'Tensor_Analysis/Pressure_Tensor_Output.dat', 'a') as fp2:
            np.savetxt(fp2, results[1], delimiter='        ', fmt='%.14F')

        #with get_lock(f'Tensor_Analysis/Gyration_Tensor_Output.dat.lock'):
        with open(f'Tensor_Analysis/Gyration_Tensor_Output.dat', 'a') as fp3:
            np.savetxt(fp3, results[2], delimiter='        ', fmt='%.14F')

        #with get_lock(f'Tensor_Analysis/Stress_Eigenvalues_Output.dat.lock'):
        with open(f'Tensor_Analysis/Stress_Eigenvalues_Output.dat', 'a') as fp4:
           np.savetxt(fp4, results[3], delimiter='        ', fmt='%.14F')

        #with get_lock(f'Tensor_Analysis/Stress_Eigenvectors_Output.dat.lock'):
        with open(f'Tensor_Analysis/Stress_Eigenvectors_Output.dat', 'a') as fp5:
            np.savetxt(fp5, results[4], delimiter='        ', fmt='%.14F')

        #with get_lock(f'Tensor_Analysis/Pressure_Eigenvalues_Output.dat.lock'):
        with open(f'Tensor_Analysis/Pressure_Eigenvalues_Output.dat', 'a') as fp6:
            np.savetxt(fp6, results[5], delimiter='        ', fmt='%.14F')
                           
        #with get_lock(f'Tensor_Analysis/Pressure_Eigenvectors_Output.dat.lock'):
        with open(f'Tensor_Analysis/Pressure_Eigenvectors_Output.dat', 'a') as fp7: 
            np.savetxt(fp7, results[6], delimiter='        ', fmt='%.14F')
        
        #with get_lock(f'Tensor_Analysis/Gyration_Eigenvalues_Output.dat.lock'):
        with open(f'Tensor_Analysis/Gyration_Eigenvalues_Output.dat', 'a') as fp8:
            np.savetxt(fp8, results[7], delimiter='        ', fmt='%.14F')
                           
        #with get_lock(f'Tensor_Analysis/Gyration_Eigenvectors_Output.dat.lock'):
        with open(f'Tensor_Analysis/Gyration_Eigenvectors_Output.dat', 'a') as fp9:
             np.savetxt(fp9, results[8], delimiter='        ', fmt='%.14F')        

        #with get_lock(f'Tensor_Analysis/Radius_Gyration_Output.dat.lock'):
        with open(f'Tensor_Analysis/Radius_Gyration_Output.dat', 'a') as fp10:
             np.savetxt(fp10, results[9], delimiter='        ', fmt='%.14F')

	#with get_lock(f'Tensor_Analysis/Angle_Max_Stress_Gyration_Output.dat.lock'):
        with open(f'Tensor_Analysis/Angle_Max_Stress_Gyration_Output.dat', 'a') as fp11:
             np.savetxt(fp11, np.array([angle]), delimiter='        ', fmt='%.14F')

        #with get_lock(f'Tensor_Analysis/Max_Eigenvalues_Output.dat.lock'):
        with open(f'Tensor_Analysis/Max_Eigenvalues_Gyration_Stress_Output.dat', 'a') as fp12:
            np.savetxt(fp12, [[max_eigenvalue_gyration, max_eigenvalue_stress]], delimiter='        ', fmt='%.14F')
        

def polymer_analysis():
    # Read the coordinates from an external file containing all time frames
    num_beads_per_polymer = 20
    coords_all = np.loadtxt("File1.dat")
    # Number of particles and frames
    num_particles = 4000
    num_frames = coords_all.shape[0] // num_particles
    print(num_frames)
    # Reshape the coordinates to create a 3D array with shape (num_time_frames, num_polymers, num_beads_per_polymer, 2)
    coords_all = coords_all.reshape(num_frames, num_particles, 2)
    num_polymers = num_particles // num_beads_per_polymer
    print(num_polymers)
    coords_all = coords_all.reshape(num_frames, num_polymers, num_beads_per_polymer, 2)
    avg_largest_eigenvalue = []
    
    num_cores = 6 #os.cpu_count() // 2
  
    args_list = [(time_frame, coords_all) for time_frame in range(num_frames)]
    job.Parallel(n_jobs=num_cores)(job.delayed(polymer_analysis_wrapper)(args) for args in args_list)

if __name__ == '__main__':
    polymer_analysis()
