# export_step.py
import uuid

from typing_extensions import Literal

from OCP.XSControl import XSControl_WorkSession
from OCP.STEPCAFControl import STEPCAFControl_Writer
from OCP.STEPControl import STEPControl_StepModelType
from OCP.IFSelect import IFSelect_ReturnStatus
from OCP.TCollection import TCollection_ExtendedString, TCollection_AsciiString
from OCP.Interface import Interface_Static
from OCP.XCAFDoc import XCAFDoc_DocumentTool

from cadquery.occ_impl.assembly import AssemblyProtocol, toCAF, toFusedCAF

class ExportModes:
    DEFAULT = "default"
    FUSED = "fused"

STEPExportModeLiterals = Literal["default", "fused"]

def export_step_ap214(
    assy: AssemblyProtocol,
    path: str,
    mode: STEPExportModeLiterals = "default",
    **kwargs
) -> bool:

    # Handle the extra settings for the STEP export
    pcurves = 1
    if "write_pcurves" in kwargs and not kwargs["write_pcurves"]:
        pcurves = 0
    precision_mode = kwargs["precision_mode"] if "precision_mode" in kwargs else 0
    fuzzy_tol = kwargs["fuzzy_tol"] if "fuzzy_tol" in kwargs else None
    glue = kwargs["glue"] if "glue" in kwargs else False

    # Use the assembly name if the user set itexportAssembly
    assembly_name = assy.name if assy.name else str(uuid.uuid1())

    # Handle the doc differently based on which mode we are using
    # creating XCAF Document and opening an Application
    if mode == "fused":
        _, doc = toFusedCAF(assy, glue, fuzzy_tol)
    else:  # Includes "default"
        _, doc = toCAF(assy, True)

    #Try to insert Notes into zhe document.
    ntool = XCAFDoc_DocumentTool.NotesTool_s(doc.Main())
    user = TCollection_ExtendedString(TCollection_AsciiString("User"))
    timestamp = TCollection_ExtendedString(TCollection_AsciiString("Timestamp"))
    comment = TCollection_ExtendedString(TCollection_AsciiString("Test text."))

    ntool.CreateComment(user, timestamp, comment)

    #STEP export
    session = XSControl_WorkSession()
    writer = STEPCAFControl_Writer(session, False)
    writer.SetColorMode(True)
    writer.SetLayerMode(True)
    writer.SetNameMode(True)
    Interface_Static.SetIVal_s("write.step.schema", 4)
    Interface_Static.SetIVal_s("write.surfacecurve.mode", pcurves)
    Interface_Static.SetIVal_s("write.precision.mode", precision_mode)
    writer.Transfer(doc, STEPControl_StepModelType.STEPControl_AsIs)
    status = writer.Write(path)

    return status == IFSelect_ReturnStatus.IFSelect_RetDone