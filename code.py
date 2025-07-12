import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def simulate_bertrands_paradox_3d(threshold_length, number_of_iterations):
    # We store the cumulative counts
    cumulative_counts_method1 = [0, 0]
    cumulative_counts_method2 = [0, 0]

    plt.ion()

    for i in range(number_of_iterations):
        plt.clf()

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        # Method 1:By randomly selecting two end points on the unit sphere
        theta1 = np.random.uniform(0, 2 * np.pi)
        phi1 = np.random.uniform(0, np.pi)
        theta2 = np.random.uniform(0, 2 * np.pi)
        phi2 = np.random.uniform(0, np.pi)
        x1, y1, z1 = np.sin(phi1) * np.cos(theta1), np.sin(phi1) * np.sin(theta1), np.cos(phi1)
        x2, y2, z2 = np.sin(phi2) * np.cos(theta2), np.sin(phi2) * np.sin(theta2), np.cos(phi2)
        chord_length_method1 = np.sqrt((x2 - x1)**2 + (y2 - y1)**2 + (z2 - z1)**2)

        # Method 2:By randomly selecting the midpoint inside the unit sphere
        r = np.random.uniform(0, 1)
        theta = np.random.uniform(0, 2 * np.pi)
        phi = np.random.uniform(0, np.pi)
        x_mid, y_mid, z_mid = r * np.sin(phi) * np.cos(theta), r * np.sin(phi) * np.sin(theta), r * np.cos(phi)

        # We should calculate the perpendicular direction using cross product
        radial_vector = np.array([x_mid, y_mid, z_mid])
        random_vector = np.random.randn(3)
        perpendicular_vector = np.cross(radial_vector, random_vector)
        perpendicular_vector = perpendicular_vector / np.linalg.norm(perpendicular_vector)

        # Now we find the length to the surface of the sphere
        chord_length_method2 = 2 * np.sqrt(1 - r**2)

        # We should calculate the endpoints of the chord
        half_length = np.sqrt(1 - r**2)
        endpoint1 = radial_vector + half_length * perpendicular_vector
        endpoint2 = radial_vector - half_length * perpendicular_vector
        x1_method2, y1_method2, z1_method2 = endpoint1
        x2_method2, y2_method2, z2_method2 = endpoint2

        # Now we need to update the cumulative counts and determine colors based on length
        if chord_length_method1 < threshold_length:
            cumulative_counts_method1[0] += 1
            color1 = 'red'
        else:
            cumulative_counts_method1[1] += 1
            color1 = 'blue'

        if chord_length_method2 < threshold_length:
            cumulative_counts_method2[0] += 1
            color2 = 'red'
        else:
            cumulative_counts_method2[1] += 1
            color2 = 'blue'

        # Now we plot the unit sphere
        u = np.linspace(0, 2 * np.pi, 100)
        v = np.linspace(0, np.pi, 100)
        x = np.outer(np.cos(u), np.sin(v))
        y = np.outer(np.sin(u), np.sin(v))
        z = np.outer(np.ones(np.size(u)), np.cos(v))
        ax.plot_surface(x, y, z, color='black', alpha=0.1)

        # We are going to plot the chords with updated colors
        ax.plot([x1, x2], [y1, y2], [z1, z2], color=color1, label="Method 1")
        ax.plot([x1_method2, x2_method2], [y1_method2, y2_method2], [z1_method2, z2_method2], color=color2, label="Method 2")

        # We should display cumulative counts and simulated probabilities
        total_counts1 = sum(cumulative_counts_method1)
        total_counts2 = sum(cumulative_counts_method2)
        ax.set_title(f"Iteration: {i + 1}\n"
                  f"Method 1: Shorter={cumulative_counts_method1[0]}, Longer={cumulative_counts_method1[1]}, "
                  f"Probability Longer={cumulative_counts_method1[1] / total_counts1:.2f}\n"
                  f"Method 2: Shorter={cumulative_counts_method2[0]}, Longer={cumulative_counts_method2[1]}, "
                  f"Probability Longer={cumulative_counts_method2[1] / total_counts2:.2f}")

        ax.set_xlabel("x-axis")
        ax.set_ylabel("y-axis")
        ax.set_zlabel("z-axis")
        ax.legend(loc="upper right")
        ax.set_xlim([-1.5, 1.5])
        ax.set_ylim([-1.5, 1.5])
        ax.set_zlim([-1.5, 1.5])
        plt.draw()
        plt.pause(0.1)

    plt.ioff()
    plt.show()

# Example run:
simulate_bertrands_paradox_3d(0.1, 100)