# Rudolf C. Kischer
# 260956107

import numpy as np
import polyscope as ps
import polyscope.imgui as psim
import argparse, json
from rigidbody import *
from contact import *
from collision import *

# Initialize polyscope, ground plane at zero
ps.init()
ps.set_automatically_compute_scene_extents(False)
ps.set_length_scale(10)
low = np.array((-2, -2., -2.)) 
high = np.array((2., 2., 2.))
ps.set_bounding_box(low, high)
ps.set_always_redraw(True)
ps.set_ground_plane_height_factor(-0.2,is_relative=True)
			
# process argumetns from command line and setup with selected json file
parser = argparse.ArgumentParser()
parser.add_argument("--file", type=str, default = "scenes/scene2.json")
args = parser.parse_args()
data = json.load(open(args.file))
rigid_body_list = []
for body_desc in data['bodies']:
	rigid_body_list.append( RigidBody( body_desc ) )
sim_params = data['sim_parameters'][0]
gravity = np.array(sim_params.get('gravity',(0.0,0.0,0.0)))	# 3x1 vector
h = sim_params.get('dt',0.01)
mu = sim_params.get('mu',0)
substeps = sim_params.get('substeps',1)
is_running = sim_params.get('is_running',False)
check_collisions = sim_params.get('check_collisions',True)
show_contacts = sim_params.get('show_contacts',False)
stop_on_contact = sim_params.get('stop_on_contact',False)

# other setup before main loop
collision = Collision()
elapsed = 0
np.set_printoptions(formatter={'float_kind':"{:.2f}".format})

def main_display_loop():
  global is_running, check_collisions, show_contacts, stop_on_contact, substeps, gravity, mu, elapsed, h

  # TODO: SET YOUR NAME AND STUDENT NUMBER HERE
  psim.TextUnformatted("Rudolf C. Kischer 260956107")
  psim.TextUnformatted(sim_params['name'])
  if(psim.Button("Reset")):
    collision.reset()
    collision.update_display(show_contacts)
    for rb in rigid_body_list:	
      rb.reset()
      rb.update_display()
      elapsed = 0
  do_step = False
  if(psim.Button("Step")):
    do_step = True
  _, is_running	= psim.Checkbox("Run", is_running) 
  psim.TextUnformatted("Elapsed = " + str(elapsed))
  _, h = psim.SliderFloat("step size", h, v_min=0.001, v_max=0.1)
  _, substeps = psim.InputInt("substeps", substeps, step=1, step_fast=1)
  if substeps < 1:
    substeps = 1; 
  _, check_collisions = psim.Checkbox("Check collisions", check_collisions )
  _, stop_on_contact = psim.Checkbox("Stop on contact", stop_on_contact)
  changed, show_contacts = psim.Checkbox("Show contacts", show_contacts) 
  if changed:
    collision.update_display(show_contacts)
  _, gravity[1] = psim.SliderFloat("gravity y", gravity[1], v_min=-10, v_max=10)
  _, mu = psim.SliderFloat("friction", mu, v_min=0, v_max=5)

  # show the energy and momentum of just the first body, only a test option
  if(psim.TreeNode("Energy and Momentum")):
    for rb in rigid_body_list:
      psim.TextUnformatted(rb.name)
      #TODO: compute and display the kinetic energy, potential energy, and linear and angular momentum of each body	
      psim.TextUnformatted("Kinetic Energy = " + str(rb.kinetic_energy()))
      psim.TextUnformatted("Potential Energy = " + str(rb.potential_energy(gravity)))
      psim.TextUnformatted("Total Energy = " + str(rb.kinetic_energy() + rb.potential_energy(gravity)))
      psim.TextUnformatted("Linear Momentum = " + str(rb.linear_momentum()))
      psim.TextUnformatted("Angular Momentum = " + str(rb.angular_momentum()))
    psim.TreePop()

  if is_running or do_step:
    # ignore substeps for now
    for i in range(substeps):
      stepsize = h/substeps
      for rb in rigid_body_list:
        rb.zero_force_and_torque()
        rb.add_force( gravity * rb.mass )
        rb.step_vel( stepsize )
      if check_collisions and collision.check(rigid_body_list) :
        collision.process(rigid_body_list,mu)
        if stop_on_contact:
          is_running = False
      for rb in rigid_body_list:
        rb.step_pos( stepsize )
      elapsed += stepsize
    for rb in rigid_body_list:

      rb.update_display()
    if show_contacts:
      collision.update_display(show_contacts)

ps.set_user_callback(main_display_loop)
ps.show()