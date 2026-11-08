# Merge-Collector-Template-Tool
A tool for creating .dxf templates for merge collectors. Either print on paper cut out and trace onto your tubing, or use with supporting tube laser cutters

Usage: py merge_calculator.py &lt;D1&gt; &lt;D2&gt; &lt;Angle&gt; &lt;N&gt; [Height]

- D1 is the diameter of the tubes you are merging into one
- D2 is the diameter of the tube the others are merging into
- Angle is the total angle (in degrees) of the merge (input angle of each runner is angle/2 from the centerline)
- N is the number of tubes your merging into 1
- [Height] is the optional section height above the cut (default 25)

Length units are arbitrary, set the units in whatever dxf viewer to the same ones you used as your input arguments for this

Example mm: py merge_calculator.py 50.8 63.5 20 2  
Example in: py merge_calculator.py 1.5 2 30 3 1
