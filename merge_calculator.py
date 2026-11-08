#Usage: py merge_calculator.py <D1> <D2> <Angle> <N> [Height]

#       D1 is the diameter of the tubes you are merging into one
#       D2 is the diameter of the tube the others are merging into
#       Angle is the total angle (in degrees) of the merge (input angle of each runner is angle/2 from the centerline)
#       N is the number of tubes your merging into 1
#       [Height] is the optional section height above the cut (default 25)

#Length units are arbitrary, set the units in whatever dxf viewer to the same ones you used as input arguments for this

#Example mm: py merge_calculator.py 50.8 63.5 20 2
#Example in: py merge_calculator.py 1.5 2 30 3 1


#I derived the length equation but couldnt figure out how to solve for theta, so I gave it to gemma 4 31b to make a script to solve it. After using onshape to model the first 2 then converting to dxf, I fed a geometric description of the cuts into gemma to alter the code to draw the dxf itself, after making some minor alterations myself this is the result. I compared the resultant dxf to one generated from onshape and they were identical except for one horizontal line which ended up being useful for alignment.

#Yes I'm fully aware that calculating the actual arc length along the ellipse is overkill and using a circular approximation would work just as well for shallower angles, but doing it this way results in near ideal geometry

#Length equation in wolfram scrypt: RIntegrate[Sqrt[Divide[Power[Sin[Φ],2],Power[Cos[α*Divide[π,180]],2]]+Power[Cos[Φ],2]],{Φ,0,2π}] - 2*RIntegrate[Sqrt[Divide[Power[Sin[Φ],2],Power[Cos[α*Divide[π,180]],2]]+Power[Cos[Φ],2]],{Φ,0,Divide[θ,2]*Divide[π,180]}]

#Length equation in latex: R \int_{0}^{2\pi} \sqrt{\frac{\sin^2(\Phi)}{\cos^2\left(\alpha \times \frac{\pi}{180}\right)} + \cos^2(\Phi)} \, d\Phi - 2R \int_{0}^{\frac{\theta}{2}\frac{\pi}{180}} \sqrt{\frac{\sin^2(\Phi)}{\cos^2\left(\alpha \times \frac{\pi}{180}\right)} + \cos^2(\Phi)} \, d\Phi

#Parametric onshape model (idk if the link will expire): https://cad.onshape.com/documents/b71db93cacd2f55f292b98aa/w/dbcb92b597f70599d3f39db4/e/e52fa541d35a83b14820e04a

import sys
import numpy as np
from scipy.integrate import quad
from scipy.optimize import minimize_scalar

def integrand(phi, alpha_deg):
    """
    Evaluates the integrand function exactly as provided in the reference.
    """
    alpha_rad = np.radians(alpha_deg)
    cos_alpha = np.cos(alpha_rad)
    
    term1 = (np.sin(phi) ** 2) / (cos_alpha ** 2)
    term2 = np.cos(phi) ** 2
    
    return np.sqrt(term1 + term2)

def solve_for_theta(D1, D2, alpha_deg, n):
    """
    Solves for theta using the exact formula logic matching your attached file.
    """
    R_s = D1 / 2
    R_l = D2 / 2
    
    L = (2 * np.pi * R_l) / n
    
    full_integral, _ = quad(integrand, 0, 2 * np.pi, args=(alpha_deg,))
    
    def objective(theta_deg):
        t_deg = float(np.asarray(theta_deg).item())
        upper_limit = np.radians(t_deg / 2)
        partial_integral, _ = quad(integrand, 0, upper_limit, args=(alpha_deg,))
        
        residual = R_s * full_integral - 2 * R_s * partial_integral - L
        return residual ** 2

    result = minimize_scalar(objective, bounds=(0.0, 360.0), method='bounded')
    
    if not result.success:
        raise ValueError("Optimization failed to converge to a valid physical angle.")
        
    return float(result.x)

def export_unwrapped_dxf(filename, D1, alpha_deg, theta_deg, n, extension_height):
    """
    Generates the DXF file representing the true perpendicular cut profile development.
    Ensures the template bounding box is a fully continuous rectangle tangent to the miter ears.
    """
    R_s = D1 / 2
    alpha_rad = np.radians(alpha_deg)
    
    theta_rad = np.radians(theta_deg)
    half_theta = theta_rad / 2
    
    # 360 degrees unwrapped resolution
    num_points = 360
    phi_vals = np.linspace(0, 2 * np.pi, num_points)
    
    raw_points = []
    
    # First pass: Calculate unshifted heights to find the true peak of the miter ears
    for phi in phi_vals:
        phi_norm = phi if phi <= np.pi else phi - 2 * np.pi
        x = R_s * phi
        
        # Base miter height profile function from your original code reference
        y_base = R_s * (1 - np.cos(half_theta)) * np.tan(alpha_rad) + R_s * (np.cos(phi_norm) - np.cos(half_theta)) / np.tan(alpha_rad)
        
        if abs(phi_norm) <= half_theta:
            if n == 2:
                y = y_base
            else:
                beta_rad = np.radians(180.0 / n)
                miter_offset = (R_s * (np.sin(half_theta) - np.sin(abs(phi_norm)))) / (np.tan(beta_rad) * np.sin(alpha_rad))
                y = y_base - miter_offset
        else:
            y = R_s * (1 - np.cos(phi_norm)) * np.tan(alpha_rad)
            
        raw_points.append((x, y))
        
    # Aligns with your exact syntax: Find the maximum height of the miter ears
    max_ears_y = max(p[1] for p in raw_points)
    
    # Shift profile points so the highest miter ear peaks are tangent exactly at Y = 0
    profile_points = [(p[0], p[1] - max_ears_y) for p in raw_points]
    
    total_width = R_s * 2 * np.pi
    box_bottom_y = 0.0
    box_top_y = extension_height
    
    with open(filename, "w") as dxf:
        dxf.write("0\nSECTION\n2\nHEADER\n0\nENDSEC\n")
        dxf.write("0\nSECTION\n2\nENTITIES\n")
        
        # 1. Write the clean, single continuous cut profile path line
        for i in range(len(profile_points) - 1):
            p1 = profile_points[i]
            p2 = profile_points[i+1]
            dxf.write("0\nLINE\n8\nCut_Profile\n")
            dxf.write(f"10\n{p1[0]:.6f}\n20\n{p1[1]:.6f}\n30\n0.0\n")
            dxf.write(f"11\n{p2[0]:.6f}\n21\n{p2[1]:.6f}\n31\n0.0\n")
            
        # FIXED: Removed the loop closing segment that previously linked profile_points[-1] to profile_points[0]
        
        # 2. Draw Bounding Rectangle Box Left Side Edge
        dxf.write("0\nLINE\n8\nTemplate_Outline\n")
        dxf.write(f"10\n0.0\n20\n{box_bottom_y:.6f}\n30\n0.0\n")
        dxf.write(f"11\n0.0\n21\n{box_top_y:.6f}\n31\n0.0\n")
        
        # 3. Draw Bounding Rectangle Box Right Side Edge
        dxf.write("0\nLINE\n8\nTemplate_Outline\n")
        dxf.write(f"10\n{total_width:.6f}\n20\n{box_bottom_y:.6f}\n30\n0.0\n")
        dxf.write(f"11\n{total_width:.6f}\n21\n{box_top_y:.6f}\n31\n0.0\n")
        
        # 4. Draw Bounding Rectangle Box Straight Flat Top Edge Line
        dxf.write("0\nLINE\n8\nTemplate_Outline\n")
        dxf.write(f"10\n0.0\n20\n{box_top_y:.6f}\n30\n0.0\n")
        dxf.write(f"11\n{total_width:.6f}\n21\n{box_top_y:.6f}\n31\n0.0\n")
        
        # 5. Draw Bounding Rectangle Box Straight Flat Bottom Edge Line (Tangent to ears at Y = 0)
        dxf.write("0\nLINE\n8\nTemplate_Outline\n")
        dxf.write(f"10\n0.0\n20\n{box_bottom_y:.6f}\n30\n0.0\n")
        dxf.write(f"11\n{total_width:.6f}\n21\n{box_bottom_y:.6f}\n31\n0.0\n")
        
        # 6. Draw Left Vertical Connector Line to close the gap up to the box (Green Line)
        dxf.write("0\nLINE\n8\nTemplate_Outline\n")
        dxf.write(f"10\n0.0\n20\n{profile_points[0][1]:.6f}\n30\n0.0\n")
        dxf.write(f"11\n0.0\n21\n{box_bottom_y:.6f}\n31\n0.0\n")
        
        # 7. Draw Right Vertical Connector Line to close the gap up to the box (Green Line)
        dxf.write("0\nLINE\n8\nTemplate_Outline\n")
        dxf.write(f"10\n{total_width:.6f}\n20\n{profile_points[-1][1]:.6f}\n30\n0.0\n")
        dxf.write(f"11\n{total_width:.6f}\n21\n{box_bottom_y:.6f}\n31\n0.0\n")
        
        dxf.write("0\nENDSEC\n0\nEOF\n")

if __name__ == "__main__":
    if len(sys.argv) < 5 or len(sys.argv) > 6:
        print("Error: Invalid number of arguments.")
        print("Usage:\tpy merge_calculator.py <D1> <D2> <Angle> <N> [Height]")
        print("\tD1 is the diameter of the tubes you are merging into one")
        print("\tD2 is the diameter of the tube the others are merging into")
        print("\tAngle is the total angle of the merge (inputs are angle/2 from the centerline)")
        print("\tN is the number of tubes your merging into 1")
        print("\t[Height] is the optional section height above the cut (default 25)")
        print("\r\nLength units are arbitraty, in whatever tool you use to view/print the dxf use the same units as you entered here")
        print("Example mm: py merge_calculator.py 50.8 63.5 20 2")
        print("Example in: py merge_calculator.py 1.5 2 30 3 1\r\n")
        sys.exit(1)
        
    try:
        D1_val    = float(sys.argv[1])
        D2_val    = float(sys.argv[2])
        alpha_val = float(sys.argv[3]) / 2  # Centerline offset angle (Angle / 2)
        n_val     = int(sys.argv[4])
        
        if len(sys.argv) == 6:
            h_raw = float(sys.argv[5])
            height_val = int(h_raw) if h_raw.is_integer() else h_raw
        else:
            height_val = 25
        
        if n_val == 0:
            print("Error: 'N' cannot be zero.")
            sys.exit(1)
            
        theta_ans = solve_for_theta(D1_val, D2_val, alpha_val, n_val)
        print(f"Solved Theta: {theta_ans:.4f} degrees")
        print(f"Using top extension height: {height_val} units")
        
        dxf_filename = f"Merge_Template_{D1_val}_{D2_val}_{float(sys.argv[3])}_{n_val}_{height_val}.dxf"
        export_unwrapped_dxf(dxf_filename, D1_val, alpha_val, theta_ans, n_val, height_val)
        print(f"Exported true geometric unwrapped DXF template to: {dxf_filename}")
        
    except ValueError as ve:
        print(f"Error: {ve}")
        sys.exit(1)
    except Exception as e:
        print(f"Execution Error: {e}")

