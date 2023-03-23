<!DOCTYPE qgis PUBLIC 'http://mrcc.com/qgis.dtd' 'SYSTEM'>
<qgis minScale="100000000" simplifyDrawingTol="1" hasScaleBasedVisibilityFlag="0" labelsEnabled="0" simplifyDrawingHints="0" version="3.16.1-Hannover" simplifyLocal="1" maxScale="0" readOnly="0" styleCategories="AllStyleCategories" simplifyMaxScale="1" simplifyAlgorithm="0">
  <flags>
    <Identifiable>1</Identifiable>
    <Removable>1</Removable>
    <Searchable>1</Searchable>
  </flags>
  <temporal enabled="0" endField="" durationField="" accumulate="0" endExpression="" fixedDuration="0" startField="" durationUnit="min" mode="0" startExpression="">
    <fixedRange>
      <start></start>
      <end></end>
    </fixedRange>
  </temporal>
  <renderer-v2 enableorderby="0" symbollevels="0" forceraster="0" type="singleSymbol">
    <symbols>
      <symbol alpha="1" clip_to_extent="1" force_rhr="0" type="marker" name="0">
        <layer class="SimpleMarker" pass="0" enabled="1" locked="0">
          <prop v="0" k="angle"/>
          <prop v="0,144,233,191" k="color"/>
          <prop v="1" k="horizontal_anchor_point"/>
          <prop v="miter" k="joinstyle"/>
          <prop v="arrow" k="name"/>
          <prop v="0,0" k="offset"/>
          <prop v="3x:0,0,0,0,0,0" k="offset_map_unit_scale"/>
          <prop v="Inch" k="offset_unit"/>
          <prop v="35,35,35,255" k="outline_color"/>
          <prop v="solid" k="outline_style"/>
          <prop v="0" k="outline_width"/>
          <prop v="3x:0,0,0,0,0,0" k="outline_width_map_unit_scale"/>
          <prop v="Inch" k="outline_width_unit"/>
          <prop v="diameter" k="scale_method"/>
          <prop v="2" k="size"/>
          <prop v="3x:0,0,0,0,0,0" k="size_map_unit_scale"/>
          <prop v="MapUnit" k="size_unit"/>
          <prop v="1" k="vertical_anchor_point"/>
          <data_defined_properties>
            <Option type="Map">
              <Option value="" type="QString" name="name"/>
              <Option type="Map" name="properties">
                <Option type="Map" name="angle">
                  <Option value="true" type="bool" name="active"/>
                  <Option value="HWD" type="QString" name="field"/>
                  <Option value="2" type="int" name="type"/>
                </Option>
                <Option type="Map" name="size">
                  <Option value="true" type="bool" name="active"/>
                  <Option value="&quot;HWS&quot;/MEDIAN(&quot;HWS&quot;)" type="QString" name="expression"/>
                  <Option value="3" type="int" name="type"/>
                </Option>
              </Option>
              <Option value="collection" type="QString" name="type"/>
            </Option>
          </data_defined_properties>
        </layer>
      </symbol>
    </symbols>
    <rotation/>
    <sizescale/>
  </renderer-v2>
  <customproperties>
    <property value="0" key="embeddedWidgets/count"/>
    <property key="variableNames"/>
    <property key="variableValues"/>
  </customproperties>
  <blendMode>0</blendMode>
  <featureBlendMode>0</featureBlendMode>
  <layerOpacity>1</layerOpacity>
  <SingleCategoryDiagramRenderer attributeLegend="1" diagramType="Histogram">
    <DiagramCategory penWidth="0" sizeScale="3x:0,0,0,0,0,0" penAlpha="255" penColor="#000000" minScaleDenominator="0" direction="0" width="15" scaleBasedVisibility="0" scaleDependency="Area" minimumSize="0" backgroundColor="#ffffff" spacingUnit="MM" labelPlacementMethod="XHeight" diagramOrientation="Up" lineSizeType="MM" lineSizeScale="3x:0,0,0,0,0,0" opacity="1" showAxis="1" maxScaleDenominator="1e+08" backgroundAlpha="255" height="15" spacing="5" barWidth="5" rotationOffset="270" spacingUnitScale="3x:0,0,0,0,0,0" sizeType="MM" enabled="0">
      <fontProperties description="Sans Serif,9,-1,5,50,0,0,0,0,0" style=""/>
      <attribute label="" color="#000000" field=""/>
      <axisSymbol>
        <symbol alpha="1" clip_to_extent="1" force_rhr="0" type="line" name="">
          <layer class="SimpleLine" pass="0" enabled="1" locked="0">
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
                <Option value="" type="QString" name="name"/>
                <Option name="properties"/>
                <Option value="collection" type="QString" name="type"/>
              </Option>
            </data_defined_properties>
          </layer>
        </symbol>
      </axisSymbol>
    </DiagramCategory>
  </SingleCategoryDiagramRenderer>
  <DiagramLayerSettings linePlacementFlags="18" zIndex="0" dist="0" placement="0" priority="0" showAll="1" obstacle="0">
    <properties>
      <Option type="Map">
        <Option value="" type="QString" name="name"/>
        <Option name="properties"/>
        <Option value="collection" type="QString" name="type"/>
      </Option>
    </properties>
  </DiagramLayerSettings>
  <geometryOptions removeDuplicateNodes="0" geometryPrecision="0">
    <activeChecks/>
    <checkConfiguration/>
  </geometryOptions>
  <legend type="default-vector"/>
  <referencedLayers/>
  <fieldConfiguration>
    <field name="ID_POINT" configurationFlags="None">
      <editWidget type="Range">
        <config>
          <Option/>
        </config>
      </editWidget>
    </field>
    <field name="HWS" configurationFlags="None">
      <editWidget type="TextEdit">
        <config>
          <Option/>
        </config>
      </editWidget>
    </field>
    <field name="HWD" configurationFlags="None">
      <editWidget type="TextEdit">
        <config>
          <Option/>
        </config>
      </editWidget>
    </field>
    <field name="VWS" configurationFlags="None">
      <editWidget type="TextEdit">
        <config>
          <Option/>
        </config>
      </editWidget>
    </field>
    <field name="WS" configurationFlags="None">
      <editWidget type="TextEdit">
        <config>
          <Option/>
        </config>
      </editWidget>
    </field>
  </fieldConfiguration>
  <aliases>
    <alias name="" field="ID_POINT" index="0"/>
    <alias name="" field="HWS" index="1"/>
    <alias name="" field="HWD" index="2"/>
    <alias name="" field="VWS" index="3"/>
    <alias name="" field="WS" index="4"/>
  </aliases>
  <defaults>
    <default expression="" applyOnUpdate="0" field="ID_POINT"/>
    <default expression="" applyOnUpdate="0" field="HWS"/>
    <default expression="" applyOnUpdate="0" field="HWD"/>
    <default expression="" applyOnUpdate="0" field="VWS"/>
    <default expression="" applyOnUpdate="0" field="WS"/>
  </defaults>
  <constraints>
    <constraint unique_strength="0" constraints="0" exp_strength="0" field="ID_POINT" notnull_strength="0"/>
    <constraint unique_strength="0" constraints="0" exp_strength="0" field="HWS" notnull_strength="0"/>
    <constraint unique_strength="0" constraints="0" exp_strength="0" field="HWD" notnull_strength="0"/>
    <constraint unique_strength="0" constraints="0" exp_strength="0" field="VWS" notnull_strength="0"/>
    <constraint unique_strength="0" constraints="0" exp_strength="0" field="WS" notnull_strength="0"/>
  </constraints>
  <constraintExpressions>
    <constraint desc="" field="ID_POINT" exp=""/>
    <constraint desc="" field="HWS" exp=""/>
    <constraint desc="" field="HWD" exp=""/>
    <constraint desc="" field="VWS" exp=""/>
    <constraint desc="" field="WS" exp=""/>
  </constraintExpressions>
  <expressionfields/>
  <attributeactions>
    <defaultAction value="{00000000-0000-0000-0000-000000000000}" key="Canvas"/>
  </attributeactions>
  <attributetableconfig sortExpression="" sortOrder="0" actionWidgetStyle="dropDown">
    <columns>
      <column width="-1" hidden="0" type="field" name="ID_POINT"/>
      <column width="-1" hidden="0" type="field" name="HWS"/>
      <column width="-1" hidden="0" type="field" name="HWD"/>
      <column width="-1" hidden="0" type="field" name="VWS"/>
      <column width="-1" hidden="0" type="field" name="WS"/>
      <column width="-1" hidden="1" type="actions"/>
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
    <field editable="1" name="HWD"/>
    <field editable="1" name="HWS"/>
    <field editable="1" name="ID_POINT"/>
    <field editable="1" name="VWS"/>
    <field editable="1" name="WS"/>
  </editable>
  <labelOnTop>
    <field name="HWD" labelOnTop="0"/>
    <field name="HWS" labelOnTop="0"/>
    <field name="ID_POINT" labelOnTop="0"/>
    <field name="VWS" labelOnTop="0"/>
    <field name="WS" labelOnTop="0"/>
  </labelOnTop>
  <dataDefinedFieldProperties/>
  <widgets/>
  <previewExpression>"ID_POINT"</previewExpression>
  <mapTip></mapTip>
  <layerGeometryType>0</layerGeometryType>
</qgis>
