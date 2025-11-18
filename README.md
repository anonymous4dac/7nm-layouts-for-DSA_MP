# 7nm-layouts-for-DSA_MP
Here’s a short GitHub-style description:  A benchmark for evaluating DSA_MP performance at the 7 nm node, featuring 15 synthetic layouts and one real layout with full test datasets and design rules provided.
## wire generation
- **wire_cd (w)**: The critical dimension of metal wire shapes, defined by wire width.  
  - Value: **21 nm**

- **track_pitch (p)**: The distance between adjacent wire tracks.  
  - Metal 1 Value: **42 nm**
  - Metal 2 Value: **31.5 nm**

- **min_t2t (t₁)**: The minimum line-end to line-end distance in each wire track.  
  - Value: **42 nm**

- **max_t2t (t₂)**: The maximum line-end to line-end distance in each wire track.  
  - Value: **315 nm**

- **t2t_grid (tg)**: The unit size of line-end to line-end distance.  
  - Value: **21 nm**

- **min_length (l₁)**: The minimum length of a single wire shape along the wire tracks.  
  - Value: **42 nm**

- **max_length (l₂)**: The maximum length of a single wire shape along the wire tracks.  
  - Value: **378 nm**

- **total_x (xₜ)**: The cell bounding box size in the x direction.  

- **total_y (yₜ)**: The cell bounding box size in the y direction.
- ## via generation

- **via_x (vₓ)**: Via size along the x direction.  
  - Value: **21 nm**

- **via_y (vᵧ)**: Via size along the y direction.  
  - Value: **21 nm**

- **density (ρ)**: The probability of a via appearing at a candidate via location.  
  - Value: **0.5**

- **enclosure_x (eₓ)**: The minimum horizontal distance from a via to a metal line-end (x direction).  
  - Value: **5 nm**

- **enclosure_y (eᵧ)**: The minimum horizontal distance from a via to a metal line-end (y direction).  
  - Value: **5 nm**

- **via_pitch_x (vpₓ)**: The minimum center-to-center distance of two vias in the same tract in the x direction.  
  - Value: **32.5 nm**

- **via_pitch_y (vpᵧ)**: The minimum center-to-center distance of two vias in the same tract in the y direction.  
  - Value: **22 nm**
