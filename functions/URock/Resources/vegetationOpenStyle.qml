<!DOCTYPE qgis PUBLIC 'http://mrcc.com/qgis.dtd' 'SYSTEM'>
<qgis maxScale="0" minScale="100000000" simplifyDrawingHints="1" version="3.16.1-Hannover" styleCategories="AllStyleCategories" simplifyLocal="1" readOnly="0" simplifyDrawingTol="1" simplifyAlgorithm="0" hasScaleBasedVisibilityFlag="0" simplifyMaxScale="1" labelsEnabled="0">
  <flags>
    <Identifiable>1</Identifiable>
    <Removable>1</Removable>
    <Searchable>1</Searchable>
  </flags>
  <temporal accumulate="0" durationField="" mode="0" endField="" enabled="0" endExpression="" startField="" durationUnit="min" fixedDuration="0" startExpression="">
    <fixedRange>
      <start></start>
      <end></end>
    </fixedRange>
  </temporal>
  <renderer-v2 type="singleSymbol" symbollevels="0" enableorderby="0" forceraster="0">
    <symbols>
      <symbol type="fill" name="0" force_rhr="0" clip_to_extent="1" alpha="0.75">
        <layer locked="0" enabled="1" class="SimpleFill" pass="0">
          <prop k="border_width_map_unit_scale" v="3x:0,0,0,0,0,0"/>
          <prop k="color" v="153,51,255,255"/>
          <prop k="joinstyle" v="bevel"/>
          <prop k="offset" v="0,0"/>
          <prop k="offset_map_unit_scale" v="3x:0,0,0,0,0,0"/>
          <prop k="offset_unit" v="MM"/>
          <prop k="outline_color" v="35,35,35,255"/>
          <prop k="outline_style" v="solid"/>
          <prop k="outline_width" v="0.26"/>
          <prop k="outline_width_unit" v="MM"/>
          <prop k="style" v="solid"/>
          <data_defined_properties>
            <Option type="Map">
              <Option type="QString" name="name" value=""/>
              <Option name="properties"/>
              <Option type="QString" name="type" value="collection"/>
            </Option>
          </data_defined_properties>
        </layer>
      </symbol>
    </symbols>
    <rotation/>
    <sizescale/>
  </renderer-v2>
  <customproperties>
    <property key="embeddedWidgets/count" value="0"/>
    <property key="variableNames"/>
    <property key="variableValues"/>
  </customproperties>
  <blendMode>0</blendMode>
  <featureBlendMode>0</featureBlendMode>
  <layerOpacity>1</layerOpacity>
  <SingleCategoryDiagramRenderer diagramType="Histogram" attributeLegend="1">
    <DiagramCategory enabled="0" direction="0" sizeType="MM" penAlpha="255" labelPlacementMethod="XHeight" width="15" minimumSize="0" scaleDependency="Area" barWidth="5" scaleBasedVisibility="0" penColor="#000000" showAxis="1" maxScaleDenominator="1e+08" lineSizeType="MM" spacingUnitScale="3x:0,0,0,0,0,0" height="15" spacingUnit="MM" diagramOrientation="Up" backgroundColor="#ffffff" spacing="5" backgroundAlpha="255" minScaleDenominator="0" sizeScale="3x:0,0,0,0,0,0" lineSizeScale="3x:0,0,0,0,0,0" rotationOffset="270" opacity="1" penWidth="0">
      <fontProperties description="Sans Serif,9,-1,5,50,0,0,0,0,0" style=""/>
      <attribute color="#000000" label="" field=""/>
      <axisSymbol>
        <symbol type="line" name="" force_rhr="0" clip_to_extent="1" alpha="1">
          <layer locked="0" enabled="1" class="SimpleLine" pass="0">
            <prop k="align_dash_pattern" v="0"/>
            <prop k="capstyle" v="square"/>
            <prop k="customdash" v="5;2"/>
            <prop k="customdash_map_unit_scale" v="3x:0,0,0,0,0,0"/>
            <prop k="customdash_unit" v="MM"/>
            <prop k="dash_pattern_offset" v="0"/>
            <prop k="dash_pattern_offset_map_unit_scale" v="3x:0,0,0,0,0,0"/>
            <prop k="dash_pattern_offset_unit" v="MM"/>
            <prop k="draw_inside_polygon" v="0"/>
            <prop k="joinstyle" v="bevel"/>
            <prop k="line_color" v="35,35,35,255"/>
            <prop k="line_style" v="solid"/>
            <prop k="line_width" v="0.26"/>
            <prop k="line_width_unit" v="MM"/>
            <prop k="offset" v="0"/>
            <prop k="offset_map_unit_scale" v="3x:0,0,0,0,0,0"/>
            <prop k="offset_unit" v="MM"/>
            <prop k="ring_filter" v="0"/>
            <prop k="tweak_dash_pattern_on_corners" v="0"/>
            <prop k="use_custom_dash" v="0"/>
            <prop k="width_map_unit_scale" v="3x:0,0,0,0,0,0"/>
            <data_defined_properties>
              <Option type="Map">
                <Option type="QString" name="name" value=""/>
                <Option name="properties"/>
                <Option type="QString" name="type" value="collection"/>
              </Option>
            </data_defined_properties>
          </layer>
        </symbol>
      </axisSymbol>
    </DiagramCategory>
  </SingleCategoryDiagramRenderer>
  <DiagramLayerSettings placement="1" dist="0" showAll="1" zIndex="0" obstacle="0" linePlacementFlags="18" priority="0">
    <properties>
      <Option type="Map">
        <Option type="QString" name="name" value=""/>
        <Option name="properties"/>
        <Option type="QString" name="type" value="collection"/>
      </Option>
    </properties>
  </DiagramLayerSettings>
  <geometryOptions removeDuplicateNodes="0" geometryPrecision="0">
    <activeChecks/>
    <checkConfiguration type="Map">
      <Option type="Map" name="QgsGeometryGapCheck">
        <Option type="double" name="allowedGapsBuffer" value="0"/>
        <Option type="bool" name="allowedGapsEnabled" value="false"/>
        <Option type="QString" name="allowedGapsLayer" value=""/>
      </Option>
    </checkConfiguration>
  </geometryOptions>
  <legend type="default-vector"/>
  <referencedLayers/>
  <fieldConfiguration>
    <field configurationFlags="None" name="ID_ZONE_VEG">
      <editWidget type="Range">
        <config>
          <Option/>
        </config>
      </editWidget>
    </field>
    <field configurationFlags="None" name="MIN_HEIGHT">
      <editWidget type="TextEdit">
        <config>
          <Option/>
        </config>
      </editWidget>
    </field>
    <field configurationFlags="None" name="MAX_HEIGHT">
      <editWidget type="TextEdit">
        <config>
          <Option/>
        </config>
      </editWidget>
    </field>
    <field configurationFlags="None" name="ATTENUATIO">
      <editWidget type="TextEdit">
        <config>
          <Option/>
        </config>
      </editWidget>
    </field>
    <field configurationFlags="None" name="ID_VEG">
      <editWidget type="Range">
        <config>
          <Option/>
        </config>
      </editWidget>
    </field>
  </fieldConfiguration>
  <aliases>
    <alias name="" index="0" field="ID_ZONE_VEG"/>
    <alias name="" index="1" field="MIN_HEIGHT"/>
    <alias name="" index="2" field="MAX_HEIGHT"/>
    <alias name="" index="3" field="ATTENUATIO"/>
    <alias name="" index="4" field="ID_VEG"/>
  </aliases>
  <defaults>
    <default applyOnUpdate="0" expression="" field="ID_ZONE_VEG"/>
    <default applyOnUpdate="0" expression="" field="MIN_HEIGHT"/>
    <default applyOnUpdate="0" expression="" field="MAX_HEIGHT"/>
    <default applyOnUpdate="0" expression="" field="ATTENUATIO"/>
    <default applyOnUpdate="0" expression="" field="ID_VEG"/>
  </defaults>
  <constraints>
    <constraint unique_strength="0" notnull_strength="0" exp_strength="0" constraints="0" field="ID_ZONE_VEG"/>
    <constraint unique_strength="0" notnull_strength="0" exp_strength="0" constraints="0" field="MIN_HEIGHT"/>
    <constraint unique_strength="0" notnull_strength="0" exp_strength="0" constraints="0" field="MAX_HEIGHT"/>
    <constraint unique_strength="0" notnull_strength="0" exp_strength="0" constraints="0" field="ATTENUATIO"/>
    <constraint unique_strength="0" notnull_strength="0" exp_strength="0" constraints="0" field="ID_VEG"/>
  </constraints>
  <constraintExpressions>
    <constraint desc="" field="ID_ZONE_VEG" exp=""/>
    <constraint desc="" field="MIN_HEIGHT" exp=""/>
    <constraint desc="" field="MAX_HEIGHT" exp=""/>
    <constraint desc="" field="ATTENUATIO" exp=""/>
    <constraint desc="" field="ID_VEG" exp=""/>
  </constraintExpressions>
  <expressionfields/>
  <attributeactions>
    <defaultAction key="Canvas" value="{00000000-0000-0000-0000-000000000000}"/>
  </attributeactions>
  <attributetableconfig actionWidgetStyle="dropDown" sortOrder="0" sortExpression="">
    <columns>
      <column type="field" name="ID_ZONE_VEG" width="-1" hidden="0"/>
      <column type="field" name="MIN_HEIGHT" width="-1" hidden="0"/>
      <column type="field" name="MAX_HEIGHT" width="-1" hidden="0"/>
      <column type="field" name="ATTENUATIO" width="-1" hidden="0"/>
      <column type="field" name="ID_VEG" width="-1" hidden="0"/>
      <column type="actions" width="-1" hidden="1"/>
    </columns>
  </attributetableconfig>
  <conditionalstyles>
    <rowstyles/>
    <fieldstyles/>
  </conditionalstyles>
  <storedexpressions/>
  <editform tolerant="1"></editform>
  <editforminit/>
  <editforminitcodesource>0</editforminitcodesource>
  <editforminitfilepath></editforminitfilepath>
  <editforminitcode><![CDATA[# -*- coding: utf-8 -*-
"""
QGIS forms can have a Python function that is called when the form is
opened.

Use this function to add extra logic to your forms.

Enter the name of the function in the "Python Init function"
field.
An example follows:
"""
from qgis.PyQt.QtWidgets import QWidget

def my_form_open(dialog, layer, feature):
	geom = feature.geometry()
	control = dialog.findChild(QWidget, "MyLineEdit")
]]></editforminitcode>
  <featformsuppress>0</featformsuppress>
  <editorlayout>generatedlayout</editorlayout>
  <editable>
    <field editable="1" name="ATTENUATIO"/>
    <field editable="1" name="ID_VEG"/>
    <field editable="1" name="ID_ZONE_VEG"/>
    <field editable="1" name="MAX_HEIGHT"/>
    <field editable="1" name="MIN_HEIGHT"/>
  </editable>
  <labelOnTop>
    <field labelOnTop="0" name="ATTENUATIO"/>
    <field labelOnTop="0" name="ID_VEG"/>
    <field labelOnTop="0" name="ID_ZONE_VEG"/>
    <field labelOnTop="0" name="MAX_HEIGHT"/>
    <field labelOnTop="0" name="MIN_HEIGHT"/>
  </labelOnTop>
  <dataDefinedFieldProperties/>
  <widgets/>
  <previewExpression>"ID_ZONE_VEG"</previewExpression>
  <mapTip></mapTip>
  <layerGeometryType>2</layerGeometryType>
</qgis>
