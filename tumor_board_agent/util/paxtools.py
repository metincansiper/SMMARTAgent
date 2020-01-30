from indra.java_vm import autoclass, cast
import heapq
from .sbgn import get_default_bbox

DEFAULT_MAX_CONVERSIONS = 20
REACTOME_NAME = 'reactome'
CTD_NAME = 'ctd'
PID_NAME = 'pid'
PANTHER_NAME = 'panther'

JHashSet = autoclass('java.util.HashSet')
JString = autoclass('java.lang.String')
JL3ToSBGNPDConverter = autoclass('org.biopax.paxtools.io.sbgn.L3ToSBGNPDConverter')
JByteArrayOutputStream = autoclass('java.io.ByteArrayOutputStream')
JStandardCharsets = autoclass('java.nio.charset.StandardCharsets')
JCompleter = autoclass('org.biopax.paxtools.controller.Completer')
JSimpleEditorMap = autoclass('org.biopax.paxtools.controller.SimpleEditorMap')
JCloner = autoclass('org.biopax.paxtools.controller.Cloner')
JBioPAXLevel = autoclass('org.biopax.paxtools.model.BioPAXLevel')
JSimpleIOHandler = autoclass('org.biopax.paxtools.io.SimpleIOHandler')
JModel = autoclass('org.biopax.paxtools.model.Model')
JConversion = autoclass('org.biopax.paxtools.model.level3.Conversion')
JByteArrayInputStream = autoclass('java.io.ByteArrayInputStream')
JBbox = autoclass('org.sbgn.bindings.Bbox')
JMarshaller = autoclass('javax.xml.bind.Marshaller')
JAXBContext = autoclass('javax.xml.bind.JAXBContext')
JBoolean = autoclass('java.lang.Boolean')

def datasource_score(jConversion):
    ds = get_datasource(jConversion)
    reverse_importance_order = [PANTHER_NAME, PID_NAME, CTD_NAME, REACTOME_NAME]

    if ds not in reverse_importance_order:
        return -1

    return reverse_importance_order.index(ds)

def convertToJSet(jIt):
    s = JHashSet();
    for e in jIt:
        s.add(e)

    return s

def get_datasource(jConversion):
    return jConversion.getDataSource().iterator().next().getName().iterator().next().lower()

def get_model_from_text(text, max_conversions=DEFAULT_MAX_CONVERSIONS):
    jHandler = JSimpleIOHandler()
    jInputStream = JByteArrayInputStream(JString(text.encode('utf-8')).getBytes())
    jModel = jHandler.convertFromOWL(jInputStream)
    jConversions = jModel.getObjects(JConversion)

    if max_conversions and jConversions.size() > max_conversions:
        jConversions = heapq.nlargest(max_conversions, list(jConversions.toArray()), key=datasource_score)
        jConversions = convertToJSet(jConversions)

        l3Factory = JBioPAXLevel.L3.getDefaultFactory()

        c = JCompleter(JSimpleEditorMap.L3)
        result = c.complete(jConversions, jModel)
        cln = JCloner(JSimpleEditorMap.L3, l3Factory)
        jModel = cln.clone(jModel, result)

    return jModel

def assign_bbox_to_glyphs(jGlyphs):
    for jGlyph in jGlyphs:
        bbox = get_default_bbox(jGlyph.getClazz())
        jBbox = JBbox()
        jBbox.setX(jBbox.x)
        jBbox.setY(jBbox.y)
        jBbox.setW(jBbox.w)
        jBbox.setH(jBbox.h)
        jGlyph.setBbox(jBbox)
        assign_bbox_to_glyphs(list(jGlyph.getGlyph().toArray()))

def bipax_model_to_sbgn(jModel):
    do_layout = False # For some reason being problematic if we set this to True
    jConverter = JL3ToSBGNPDConverter(None, None, do_layout)
    jBaos = JByteArrayOutputStream()
    output_stream = cast('java.io.OutputStream', jBaos)

    # jConverter.writeSBGN method does not create the bbox objects for
    # the glyphs, not sure but maybe it would be creating them if we were able
    # to run it with layout enabled. Therefore, we need to jConverter.createSBGN method
    # and do the rest of the work below.
    jSbgn = jConverter.createSBGN(jModel)
    jMap = jSbgn.getMap()
    jGlyphs = jMap.getGlyph()
    jGlyphs = list(jGlyphs.toArray())

    assign_bbox_to_glyphs(jGlyphs)

    context = JAXBContext.newInstance('org.sbgn.bindings')
    marshaller = context.createMarshaller()
    marshaller.setProperty(JMarshaller.JAXB_FORMATTED_OUTPUT, JBoolean(True))
    marshaller.marshal(jSbgn, output_stream)

    res = jBaos.toString(JStandardCharsets.UTF_8.name())

    return res

def biopax_text_to_sbgn(text, max_conversions=DEFAULT_MAX_CONVERSIONS):
    jModel = get_model_from_text(text, max_conversions)
    sbgn_text = bipax_model_to_sbgn(jModel)
    return sbgn_text
