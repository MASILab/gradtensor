import numpy as np
import scipy.io as sio
from dipy.core.sphere import disperse_charges, Sphere, HemiSphere

#n_pts = 100
theta = np.linspace(0,np.pi ,100)
sio.savemat('theta_sphere100.mat',{'theta':theta})

phi =  np.linspace(0,2 * np.pi,100)
sio.savemat('phi_sphere100.mat',{'phi':phi})

t_100 = []
p_100 = []
for t in theta:
    for p in phi:
        t_100.append(t)
        p_100.append(p)
        
#for phi in phi_all:
#    for theta in theta_all:
hsph_initial = HemiSphere(theta=t_100, phi=p_100)
hsph_updated, potential = disperse_charges(hsph_initial, 5000)

from dipy.viz import window, actor

# Enables/disables interactive visualization
interactive = False

scene = window.Scene()
scene.SetBackground(1, 1, 1)
scene.add(actor.point(hsph_initial.vertices, window.colors.red,
                      point_radius=0.05))
scene.add(actor.point(hsph_updated.vertices, window.colors.green,
                      point_radius=0.05))

print('Saving illustration as initial_vs_updated.png')
window.record(scene, out_path='initial_vs_updated.png', size=(300, 300))
if interactive:
    window.show(scene)

sph = Sphere(xyz=np.vstack((hsph_updated.vertices, -hsph_updated.vertices)))
    #print(sph)
vertices = sph.vertices
sio.savemat('vertices_sphere100.mat',{'vertices':vertices})
scene.clear()
scene.add(actor.point(sph.vertices, window.colors.green, point_radius=0.05))
    
print('Saving illustration as full_sphere.png')
window.record(scene, out_path='full_sphere.png', size=(300, 300))
if interactive:
    window.show(scene)

#from dipy.core.geometry import cart2sphere, sphere2cart, vector_norm
#for i in range(554):
#    print(cart2sphere(vertices[i,0],vertices[i,1],vertices[i,2]))
