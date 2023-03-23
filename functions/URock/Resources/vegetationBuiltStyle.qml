<!DOCTYPE qgis PUBLIC 'http://mrcc.com/qgis.dtd' 'SYSTEM'>
<qgis readOnly="0" maxScale="0" hasScaleBasedVisibilityFlag="0" labelsEnabled="0" simplifyDrawingHints="1" minScale="100000000" simplifyMaxScale="1" simplifyDrawingTol="1" version="3.16.1-Hannover" styleCategories="AllStyleCategories" simplifyAlgorithm="0" simplifyLocal="1">
  <flags>
    <Identifiable>1</Identifiable>
    <Removable>1</Removable>
    <Searchable>1</Searchable>
  </flags>
  <temporal endField="" startExpression="" accumulate="0" durationUnit="min" enabled="0" durationField="" endExpression="" mode="0" startField="" fixedDuration="0">
    <fixedRange>
      <start></start>
      <end></end>
    </fixedRange>
  </temporal>
  <renderer-v2 type="singleSymbol" forceraster="0" enableorderby="0" symbollevels="0">
    <symbols>
      <symbol clip_to_extent="1" name="0" type="fill" force_rhr="0" alpha="0.75">
        <layer locked="0" class="SimpleFill" enabled="1" pass="0">
          <prop v="3x:0,0,0,0,0,0" k="border_width_map_unit_scale"/>
          <prop v="153,153,255,255" k="color"/>
          <prop v="bevel" k="joinstyle"/>
          <prop v="0,0" k="offset"/>
          <prop v="3x:0,0,0,0,0,0" k="offset_map_unit_scale"/>
          <prop v="MM" k="offset_unit"/>
          <prop v="35,35,35,255" k="outline_color"/>
          <prop v="solid" k="outline_style"/>
          <prop v="0.26" k="outline_width"/>
          <prop v="MM" k="outline_width_unit"/>
          <prop v="solid" k="style"/>
          <data_defined_properties>
            <Option type="Map">
              <Option name="name" type="QString" value=""/>
              <Option name="properties"/>
              <Option name="type" type="QString" value="collection"/>
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
    <DiagramCategory lineSizeScale="3x:0,0,0,0,0,0" opacity="1" spacingUnitScale="3x:0,0,0,0,0,0" diagramOrientation="Up" rotationOffset="270" minScaleDenominator="0" penAlpha="255" maxScaleDenominator="1e+08" penColor="#000000" labelPlacementMethod="XHeight" direction="0" lineSizeType="MM" scaleDependency="Area" spacing="5" backgroundColor="#ffffff" penWidth="0" sizeScale="3x:0,0,0,0,0,0" enabled="0" backgroundAlpha="255" minimumSize="0" scaleBasedVisibility="0" height="15" spacingUnit="MM" showAxis="1" sizeType="MM" barWidth="5" width="15">
      <fontProperties description="Sans Serif,9,-1,5,50,0,0,0,0,0" style=""/>
      <axisSymbol>
        <symbol clip_to_extent="1" name="" type="line" force_rhr="0" alpha="1">
          <layer locked="0" class="SimpleLine" enabled="1" pass="0">
            <prop v="0" k="align_dash_pattern"/>
            <prop v="square" k="capstyle"/>
            <prop v="5;2" k="customdash"/>
            <prop v="3x:0,0,0,0,0,0" k="customdash_map_unit_scale"/>
            <prop v="MM" k="customdash_unit"/>
            <prop v="0" k="dash_pattern_offset"/>
            <prop v="3x:0,0,0,0,0,0" k="dash_pattern_offset_map_unit_scale"/>
            <prop v="MM" k="dash_pattern_offset_unit"/>
            <prop v="0" k="draw_inside_polygon"/>
            <prop v="bevel" k="joinstyle"/>
            <prop v="35,35,35,255" k="line_color"/>
            <prop v="solid" k="line_style"/>
            <prop v="0.26" k="line_width"/>
            <prop v="MM" k="line_width_unit"/>
            <prop v="0" k="offset"/>
            <prop v="3x:0,0,0,0,0,0" k="offset_map_unit_scale"/>
            <prop v="MM" k="offset_unit"/>
            <prop v="0" k="ring_filter"/>
            <prop v="0" k="tweak_dash_pattern_on_corners"/>
            <prop v="0" k="use_custom_dash"/>
            <prop v="3x:0,0,0,0,0,0" k="width_map_unit_scale"/>
            <data_defined_properties>
              <Option type="Map">
                <Option name="name" type="QString" value=""/>
                <Option name="properties"/>
                <Option name="type" type="QString" value="collection"/>
              </Option>
            </data_defined_properties>
          </layer>
        </symbol>
      </axisSymbol>
    </DiagramCategory>
  </SingleCategoryDiagramRenderer>
  <DiagramLayerSettings dist="0" linePlacementFlags="18" placement="1" showAll="1" priority="0" obstacle="0" zIndex="0">
    <properties>
      <Option type="Map">
        <Option name="name" type="QString" value=""/>
        <Option name="properties"/>
        <Option name="type" type="QString" value="collection"/>
      </Option>
    </properties>
  </DiagramLayerSettings>
  <geometryOptions geometryPrecision="0" removeDuplicateNodes="0">
    <activeChecks/>
    <checkConfiguration type="Map">
      <Option name="QgsGeometryGapCheck" type="Map">
        <Option name="allowedGapsBuffer" type="double" value="0"/>
        <Option name="allowedGapsEnabled" type="bool" value="false"/>
        <Option name="allowedGapsLayer" type="QString" value=""/>
      </Option>
    </checkConfiguration>
  </geometryOptions>
  <legend type="default-vector"/>
  <referencedLayers/>
  <fieldConfiguration>
    <field name="ID_ZONE_VEG" configurationFlags="None">
      <editWidget type="Range">
        <config>
          <Option/>
        </config>
      </editWidget>
    </field>
    <field name="MIN_HEIGHT" configurationFlags="None">
      <editWidget type="TextEdit">
        <config>
          <Option/>
        </config>
      </editWidget>
    </field>
    <field name="MAX_HEIGHT" configurationFlags="None">
      <editWidget type="TextEdit">
        <config>
          <Option/>
        </config>
      </editWidget>
    </field>
    <field name="ATTENUATIO" configurationFlags="None">
      <editWidget type="TextEdit">
        <config>
          <Option/>
        </config>
      </editWidget>
    </field>
    <field name="ID_VEG" configurationFlags="None">
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
    <default expression="" applyOnUpdate="0" field="ID_ZONE_VEG"/>
    <default expression="" applyOnUpdate="0" field="MIN_HEIGHT"/>
    <default expression="" applyOnUpdate="0" field="MAX_HEIGHT"/>
    <default expression="" applyOnUpdate="0" field="ATTENUATIO"/>
    <default expression="" applyOnUpdate="0" field="ID_VEG"/>
  </defaults>
  <constraints>
    <constraint notnull_strength="0" unique_strength="0" field="ID_ZONE_VEG" constraints="0" exp_strength="0"/>
    <constraint notnull_strength="0" unique_strength="0" field="MIN_HEIGHT" constraints="0" exp_strength="0"/>
    <constraint notnull_strength="0" unique_strength="0" field="MAX_HEIGHT" constraints="0" exp_strength="0"/>
    <constraint notnull_strength="0" unique_strength="0" field="ATTENUATIO" constraints="0" exp_strength="0"/>
    <constraint notnull_strength="0" unique_strength="0" field="ID_VEG" constraints="0" exp_strength="0"/>
  </constraints>
  <constraintExpressions>
    <constraint exp="" desc="" field="ID_ZONE_VEG"/>
    <constraint exp="" desc="" field="MIN_HEIGHT"/>
    <constraint exp="" desc="" field="MAX_HEIGHT"/>
    <constraint exp="" desc="" field="ATTENUATIO"/>
    <constraint exp="" desc="" field="ID_VEG"/>
  </constraintExpressions>
  <expressionfields/>
  <attributeactions>
    <defaultAction key="Canvas" value="{00000000-0000-0000-0000-000000000000}"/>
  </attributeactions>
  <attributetableconfig sortOrder="0" sortExpression="" actionWidgetStyle="dropDown">
    <columns>
      <column width="-1" name="ID_ZONE_VEG" type="field" hidden="0"/>
      <column width="-1" name="MIN_HEIGHT" type="field" hidden="0"/>
      <column width="-1" name="MAX_HEIGHT" type="field" hidden="0"/>
      <column width="-1" name="ATTENUATIO" type="field" hidden="0"/>
      <column width="-1" name="ID_VEG" type="field" hidden="0"/>
      <column width="-1" type="actions" hidden="1"/>
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
