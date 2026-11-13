"""STEP files (for SolidWorks, CATIA, ANSYS) built with the OpenCascade CAD kernel.

    write_step(path, xy_mm)                 2D: the closed airfoil profile as a flat face (z = 0)
    write_step(path, xy_mm, span_mm=203.2)  3D: that face extruded along z, centred on z = 0
    check_step(path)                        read a STEP file back and measure it

xy_mm is one point loop in millimetres, in Selig order (TE -> upper -> LE ->
lower -> TE), e.g. surfaces.loop(S) * chord_mm. The profile is a smooth
curve through every point (like SolidWorks' "Curve Through XYZ Points"); a
blunt trailing edge is closed with a straight line.

Needs pythonocc-core (conda-forge), listed in environment.yml. The rest of
morphgeom works without it.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
from OCC.Core.BRep import BRep_Tool
from OCC.Core.BRepBndLib import brepbndlib
from OCC.Core.BRepBuilderAPI import (BRepBuilderAPI_MakeEdge, BRepBuilderAPI_MakeFace,
                                     BRepBuilderAPI_MakeWire, BRepBuilderAPI_Transform)
from OCC.Core.BRepCheck import BRepCheck_Analyzer
from OCC.Core.BRepGProp import brepgprop
from OCC.Core.BRepMesh import BRepMesh_IncrementalMesh
from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakePrism
from OCC.Core.Bnd import Bnd_Box
from OCC.Core.GeomAPI import GeomAPI_Interpolate
from OCC.Core.GProp import GProp_GProps
from OCC.Core.IFSelect import IFSelect_RetDone
from OCC.Core.Interface import Interface_Static
from OCC.Core.Message import Message_Fail, message
from OCC.Core.STEPCAFControl import STEPCAFControl_Writer
from OCC.Core.STEPControl import STEPControl_AsIs, STEPControl_Reader
from OCC.Core.TCollection import TCollection_ExtendedString
from OCC.Core.TColgp import TColgp_HArray1OfPnt
from OCC.Core.TDataStd import TDataStd_Name
from OCC.Core.TDocStd import TDocStd_Document
from OCC.Core.TopAbs import TopAbs_FACE, TopAbs_REVERSED, TopAbs_SOLID
from OCC.Core.TopExp import TopExp_Explorer
from OCC.Core.TopLoc import TopLoc_Location
from OCC.Core.TopoDS import topods
from OCC.Core.XCAFDoc import XCAFDoc_DocumentTool
from OCC.Core.gp import gp_Pnt, gp_Trsf, gp_Vec

for _printer in message.DefaultMessenger().Printers():   # print OpenCascade errors only, not its progress chatter
    _printer.SetTraceLevel(Message_Fail)


def profile_face(xy_mm: np.ndarray):
    """The closed profile as a flat face in the XY plane."""
    xy = np.asarray(xy_mm, float)
    closed = np.allclose(xy[0], xy[-1], atol=1e-9)
    pts = xy[:-1] if closed else xy           # a closed TE: the end point is the start point
    arr = TColgp_HArray1OfPnt(1, pts.shape[0] + (1 if closed else 0))
    for i, (x, y) in enumerate(pts, start=1):
        arr.SetValue(i, gp_Pnt(float(x), float(y), 0.0))
    if closed:
        arr.SetValue(pts.shape[0] + 1, gp_Pnt(float(pts[0, 0]), float(pts[0, 1]), 0.0))
    interp = GeomAPI_Interpolate(arr, False, 1e-7)   # smooth curve through every point
    interp.Perform()
    if not interp.IsDone():
        raise RuntimeError("Could not fit a curve through the profile points.")
    wire = BRepBuilderAPI_MakeWire()
    wire.Add(BRepBuilderAPI_MakeEdge(interp.Curve()).Edge())
    if not closed:                             # blunt TE: straight line from the lower to the upper corner
        p_lo = gp_Pnt(float(xy[-1, 0]), float(xy[-1, 1]), 0.0)
        p_up = gp_Pnt(float(xy[0, 0]), float(xy[0, 1]), 0.0)
        wire.Add(BRepBuilderAPI_MakeEdge(p_lo, p_up).Edge())
    if not wire.IsDone():
        raise RuntimeError("The profile does not close.")
    face = BRepBuilderAPI_MakeFace(wire.Wire(), True)
    if not face.IsDone():
        raise RuntimeError("Could not make a flat face from the profile.")
    return face.Face()


def write_step(path, xy_mm: np.ndarray, span_mm: float = 0.0, name: str = "airfoil") -> None:
    """Write the profile as a STEP file (AP214, millimetres).

    span_mm = 0: a flat face (surface body) at z = 0, for 2D CFD.
    span_mm > 0: a solid, the face extruded along z from -span/2 to +span/2.
    """
    shape = profile_face(xy_mm)
    if span_mm > 0:
        move = gp_Trsf()
        move.SetTranslation(gp_Vec(0.0, 0.0, -span_mm / 2))
        shape = BRepBuilderAPI_Transform(shape, move, True).Shape()
        shape = BRepPrimAPI_MakePrism(shape, gp_Vec(0.0, 0.0, float(span_mm))).Shape()
    # Written through a named document so CAD shows the part as `name`
    # (the plain writer appends a counter: "name 2", "name 3", ...).
    doc = TDocStd_Document("step")
    label = XCAFDoc_DocumentTool.ShapeTool(doc.Main()).AddShape(shape, False)
    TDataStd_Name.Set(label, TCollection_ExtendedString(name))
    Interface_Static.SetCVal("write.step.schema", "AP214IS")
    Interface_Static.SetCVal("write.step.unit", "MM")
    writer = STEPCAFControl_Writer()
    writer.SetNameMode(True)
    writer.Transfer(doc, STEPControl_AsIs)
    if writer.Write(str(Path(path))) != IFSelect_RetDone:
        raise RuntimeError(f"Could not write {path}.")


def check_step(path) -> dict:
    """Read a STEP file back and measure it: is the geometry valid, its size,
    the area of each face (precise integration), and for a solid its volume.

    The volume comes from meshing the solid into small triangles: OpenCascade's
    direct volume integration is unreliable on a curve with hundreds of
    control points (it varies by tens of percent with its precision setting).
    """
    reader = STEPControl_Reader()
    if reader.ReadFile(str(Path(path))) != IFSelect_RetDone:
        raise RuntimeError(f"Could not read {path}.")
    reader.TransferRoots()
    shape = reader.OneShape()
    box = Bnd_Box()
    box.SetGap(0.0)
    brepbndlib.AddOptimal(shape, box, False, False)       # tight box on the curves themselves
    xmin, ymin, zmin, xmax, ymax, zmax = box.Get()

    faces = []
    exp = TopExp_Explorer(shape, TopAbs_FACE)
    while exp.More():
        faces.append(topods.Face(exp.Current()))
        exp.Next()
    areas = []
    for f in faces:
        props = GProp_GProps()
        brepgprop.SurfaceProperties(f, props, 1e-9)
        areas.append(props.Mass())
    n_solids = 0
    exp = TopExp_Explorer(shape, TopAbs_SOLID)
    while exp.More():
        n_solids += 1
        exp.Next()
    out = {"valid": BRepCheck_Analyzer(shape).IsValid(), "solids": n_solids, "faces": len(faces),
           "bbox_mm": (xmin, ymin, zmin, xmax, ymax, zmax), "face_areas_mm2": areas}
    if n_solids:
        out["volume_mm3"] = _mesh_volume(shape, faces)
    return out


def _mesh_volume(shape, faces) -> float:
    """Volume of a closed solid from a fine triangle mesh of its faces."""
    BRepMesh_IncrementalMesh(shape, 0.002, False, 0.05, True)
    vol = 0.0
    for f in faces:
        loc = TopLoc_Location()
        tri = BRep_Tool.Triangulation(f, loc)
        trsf = loc.Transformation()
        flip = f.Orientation() == TopAbs_REVERSED
        for i in range(1, tri.NbTriangles() + 1):
            a, b, c = tri.Triangle(i).Get()
            if flip:
                b, c = c, b
            p = [np.array(tri.Node(k).Transformed(trsf).Coord()) for k in (a, b, c)]
            vol += np.dot(p[0], np.cross(p[1], p[2])) / 6
    return abs(vol)
