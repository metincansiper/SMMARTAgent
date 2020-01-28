from indra.java_vm import autoclass, cast
import heapq

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

def datasource_score(jConversion):
    ds = get_datasource(jConversion)
    reverse_importance_order = [PANTHER_NAME, PID_NAME, CTD_NAME, REACTOME_NAME]

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

def bipax_model_to_sbgn(jModel):
    do_layout = False # For some reason being problematic if we set this to True
    jConverter = JL3ToSBGNPDConverter(None, None, do_layout)
    jBaos = JByteArrayOutputStream()
    output_stream = cast('java.io.OutputStream', jBaos)
    jConverter.writeSBGN(jModel, output_stream)
    res = jBaos.toString(JStandardCharsets.UTF_8.name())

    return res

def biopax_text_to_sbgn(text, max_conversions=DEFAULT_MAX_CONVERSIONS):
    jModel = get_model_from_text(text, max_conversions)
    sbgn_text = bipax_model_to_sbgn(jModel)
    return sbgn_text
