import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist
from sklearn.cluster import DBSCAN
import matplotlib.pyplot as plt
from scipy.spatial.distance import squareform, pdist
import multiprocessing
import os
import time

def process_frame(frame_data):
    frame_num, positions = frame_data
    
    
    # Periodic boundary conditions
    positions = np.where(positions < -35, positions + 70, positions)
    positions = np.where(positions > 35, positions - 70, positions)

    # Calculate pairwise distances between particles
    distances = squareform(pdist(positions))

    # Define the cluster cutoff
    max_distance = 1.122

    # Determine clusters using DBSCAN algorithm
    db = DBSCAN(eps=max_distance, min_samples=4, metric='precomputed').fit(distances)
    labels = db.labels_
    
    # Remove clusters with fewer than 5 particles
    unique_labels, label_counts = np.unique(labels, return_counts=True)
    for label, count in zip(unique_labels, label_counts):
    	if label != -1 and count < 4:
        	labels[labels == label] = -1

    # Count number and size of clusters
    n_clusters = len(set(labels)) - (1 if -1 in labels else 0)  
    cluster_sizes = [np.sum(labels == i) for i in range(n_clusters)]  
    Size_clusters = pd.Series([size for size in cluster_sizes if size > 0])
    n_clusters = len(Size_clusters) 
    pd.set_option('display.max_rows', None)
    
    # Create a dictionary of cluster sizes for each label
    size_dict = {}
    for label, size in zip(labels, Size_clusters):
    	if label not in size_dict:
    	    size_dict[label] = size
    	    
   # Write size_dict to a file
   #  with open(f'Clustering/Cluster_F20_S7_SizeDict_Frame{frame_num}.txt', 'w') as f:
   #  	for label, size in size_dict.items():
   #  	    f.write(f"{label}: {size}\n")
    	    
    # Remove particles that are not assigned to any cluster
    positions_clustered = positions[labels != -1]
    labels_clustered = labels[labels != -1]
    
    
    # Remove clusters with fewer than 5 particles
    unique_labels, label_counts = np.unique(labels_clustered, return_counts=True)
    for label, count in zip(unique_labels, label_counts):
    	if label != -1 and count < 4:
        	labels_clustered[labels_clustered == label] = -1
        	
    Size_clusters_pl = Size_clusters.reindex(labels_clustered)

    fig, ax = plt.subplots(figsize=(6.4,4.8), dpi=200)
    plt.subplots_adjust(left=0.1, right=0.9, bottom=0.1, top=0.9)
    plt.xlim([-35,35]), plt.ylim([-35,35])
    plt.xticks(range(-35, 36, 10))
    plt.yticks(range(-35, 36, 10))
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.rcParams.update({'font.size': 16})
    scatter = ax.scatter(positions_clustered[:,0], positions_clustered[:,1], c= Size_clusters_pl.values, alpha=1, cmap='plasma_r', marker='*', s=80, edgecolors='black', linewidths=0.05, vmin=4, vmax=max(Size_clusters))
    ax.set_xlabel('X', fontsize=16) 
    ax.set_ylabel('Y', fontsize=16)
    #ax.set_title(f'Particle Clustering_F20_Frame{frame_num}')
    cbar = fig.colorbar(scatter, shrink=0.95)
    cbar.set_label(r'$S_{clusters}$', fontsize=16, labelpad=-40, y =1.1, rotation='horizontal')
    plt.tight_layout()
    plt.savefig(f'Clustering/Clustering_F20_3F_S7_Frame{frame_num}.png')
    plt.close()

    # Print the size of each cluster and number of clusters
    with get_lock(f'Clustering/Cluster_Output_Frame{frame_num}.dat.lock'):
        with open(f'Clustering/Cluster_Output_Frame{frame_num}.dat', 'w') as fp1:
            print(f"{n_clusters}", file=fp1)

    with get_lock(f'Clustering/Cluster_Output_Frame_Frame{frame_num}.dat.lock'):
        with open(f'Clustering/Cluster_Output_Frame_Frame{frame_num}.dat', 'w') as fp2:
            print(f"{Size_clusters}", file=fp2)

    with open("Clustering/progress_S7.txt", "a") as fp3:
    	print(f"Frame {frame_num} done.", file=fp3)
    
def get_lock(lock_file, max_wait=60):
    """Acquire a lock by creating a lock file, or wait until the lock is released."""
    lock_file = lock_file.strip('.lock') + '.lock'
    start_time = time.time()
    while os.path.exists(lock_file):
        if time.time() - start_time > max_wait:
            raise TimeoutError(f"Failed to acquire lock on {lock_file}")
        time.sleep(0.1)
    open(lock_file, 'w').close()  # create lock file
    return Lock(lock_file)


class Lock:
    """A lock that releases the lock file on exit."""
    def __init__(self, lock_file):
        self.lock_file = lock_file

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        os.remove(self.lock_file)    

if __name__ == '__main__':
    #coordinates data
    coordinates = np.loadtxt("File1.dat")

    # Parallel Run 
    num_processes = 6 
    pool = multiprocessing.Pool(num_processes)

    # Number of particles and frames
    num_particles = 4000
    num_frames = coordinates.shape[0] // num_particles

    # Split coordinates into frames
    frames = [(i, coordinates[i*num_particles:(i+1)*num_particles, :]) for i in range(num_frames)]

    results = pool.map_async(process_frame, frames)

    while not results.ready():
        # do something else while waiting for the results
       # print("Waiting for results...")

    	processed_frames = results.get()

    # Close multiprocessing pool
    pool.close()
    pool.join()

    print("Done.")
