# Rotors Verification Notes

I realized that I probably need to have notes on the verification process where I log what's going on. 



## Initial Conditions
### 7/10/23
I noticed that my initial conditions were off in my aerodynamic solution. The aero solution was starting at the steady BEM solution, then sidestepping into the dynamic stall solution. This was creating a slight offset in the initial solution, which problematic for the entire solution. 

I changed to getting the intial dynamic stall solution at the initial environmental conditions, then starting from that point for the solution. So I run the BEM to get an inflow angle at t0, then I run the first time step from the starting DS states to get my initial DS states. It seems a little funky, but it more closely matches what I expect and more closely matches what OpenFAST gets. 



## Angular Deflections
### 7/11/23
Today I asked myself again if OpenFAST couples the angular deflections into the angle of attack. I remember having looked into this before, and I couldn't remember if I had found a satisfactory answer or not. 

So first off, the AeroDyn v15 manual has some information that looks important: 

<blockquote>
When AeroDyn is coupled to FAST, AeroDyn receives the instantaneous (possibly displaced/deflected) structural position, orientation, and velocities of analysis nodes in the tower, hub, and blades. As with curvature and sweep, the 2D cross sections where the blade aerodynamic analysis takes place will follow the out-of-plane deflection, but in-plane deflection is assumed to be accomplished by shearing, rather than rotation of the 2D cross section. AeroDyn also receives the local freestream (undisturbed) fluid velocities at the tower and blade nodes. (Fluid and structural calculations take place outside of the AeroDyn module and are passed as inputs to AeroDyn by the driver code.) The fluid and structural motions are provided at each coupling time step and then AeroDyn computes the aerodynamic loads on the blade and tower nodes and returns them back to FAST as part of the aero-elastic calculation. In standalone mode, the inputs to AeroDyn are prescribed by a simple driver code, without aero-elastic coupling.
</blockquote>

From this it is a little confusing on what rotations are being included and which ones are being ignored. Like, I feel like I should be able to just understand that... but I don't. I have to really think about it. The way that I read it is that flatwise deflections (deflections in my X direction, which are deflections that cause the blade to cone) are accounted for in the rotation (both linear and angular), but edgewise deflections (deflections in the structural Y direction, which are deflections that cause the blade to sweep) are only accounted for linearly (because the are modeled by a shearing action). But this says nothing about if they are rotated about radial axis. 

Previously while talking to Dr. Ning, he didn't believe me that it doesn't have rotational coupling. Probably one because he knows better, but his argument was that it would make the model useless because the angular coupling is a large portion of the dynamic. 

While reading on the NREL forum I found the following, which demonstrates that there is angular coupling that affects the angle of attack. 

<blockquote>
Jason.Jonkman - Mar '17

Dear Daniel,

The rotations sent from BeamDyn to AeroDyn, and the associated angle of attack, are correct. The rotations computed within the BeamDyn equations of motion are also correct. The error is in the section of code related to how the internal rotations are converted before being written to the output file.

Best regards,
</blockquote>

So evidently there is a coupling, and I'm pretty sure it is accomplished by rotating the loads. I've already been including the angle of attack rotations... so nothing to change there. And I'm fairly confident that I'm doing it correctly. 

So I added in the rotations due to sweep. It didn't really change anything... so yeah... I'm going to leave it for now. Hopefully those rotations were done correctly, so... yeah. 


## Divergence of loads
I was just analyzing the outputs of a 1 second solution, and I noticed the following: 
- First there is a small divergence in the loads just past 0.5 s. The model appears to recover from this, but I'm going to guess it was a difference in the dynamic stall model, considering that it appears to be a near discontinuous change (in both OF and Rotors). 
- This slight difference in the DS states causes a divergence in the inflow angle solution, which affects the DS model further (across time), which results in an overall difference in loads, which results in a difference in the deflection and other structural states, and boom... now we're at different states and there ain't now way that we're ever going to be the same at this point. 

I don't really know if this helps me much in the long run, but... it's where I'm at. 