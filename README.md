# CaliceEventDisplay
Event display program allows you to see the structure of observed hit points in 3D.

![Hit Map](img/hitmap.png?raw=true "Title")

In this program, you can
 - Read through each event reconstructed by the event building macro. ([SiWECAL-TB-monitoring](https://github.com/SiWECAL-TestBeam/SiWECAL-TB-monitoring))
   - One can have a handle on event by event analysis.
   - Currently make coincidence of `nhit_slab >= 13`.
 - Access each hit information
   - Hover curser over the hit marker. This gives you information on those hits.
   - Currently returns `hit_adc_high`, `hit_energy`, `hit_isHit`, and (`hit_x`,`hit_y`,`hit_z`,`hit_sca`)

     ![Hit Info](img/hitinfo.png?raw=true "Title")
     
 - Turn on/off the detector geometry overlay
   - One can turn off the detector geometry from Eve tab and uncheck `Geometry Scene`
   
     ![No Geometry](img/no_geometry.png?raw=true "Title")
