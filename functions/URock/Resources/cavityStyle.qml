<!DOCTYPE qgis PUBLIC 'http://mrcc.com/qgis.dtd' 'SYSTEM'>
<qgis maxScale="0" labelsEnabled="0" version="3.16.1-Hannover" simplifyDrawingTol="1" styleCategories="AllStyleCategories" hasScaleBasedVisibilityFlag="0" simplifyLocal="1" simplifyMaxScale="1" readOnly="0" minScale="100000000" simplifyDrawingHints="1" simplifyAlgorithm="0">
  <flags>
    <Identifiable>1</Identifiable>
    <Removable>1</Removable>
    <Searchable>1</Searchable>
  </flags>
  <temporal enabled="0" endField="" durationField="" startExpression="" endExpression="" durationUnit="min" fixedDuration="0" accumulate="0" mode="0" startField="">
    <fixedRange>
      <start></start>
      <end></end>
    </fixedRange>
  </temporal>
  <renderer-v2 type="singleSymbol" forceraster="0" symbollevels="0" enableorderby="0">
    <symbols>
      <symbol clip_to_extent="1" type="fill" force_rhr="0" name="0" alpha="0.75">
        <layer enabled="1" class="SimpleFill" pass="0" locked="0">
          <prop v="3x:0,0,0,0,0,0" k="border_width_map_unit_scale"/>
          <prop v="204,255,0,255" k="color"/>
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
              <Option value="" type="QString" name="name"/>
              <Option name="properties"/>
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
  <SingleCategoryDiagramRenderer diagramType="Histogram" attributeLegend="1">
    <DiagramCategory direction="0" penWidth="0" scaleDependency="Area" scaleBasedVisibility="0" spacingUnit="MM" enabled="0" spacingUnitScale="3x:0,0,0,0,0,0" sizeScale="3x:0,0,0,0,0,0" spacing="5" lineSizeScale="3x:0,0,0,0,0,0" rotationOffset="270" minScaleDenominator="0" sizeType="MM" lineSizeType="MM" backgroundColor="#ffffff" backgroundAlpha="255" penColor="#000000" maxScaleDenominator="1e+08" diagramOrientation="Up" width="15" height="15" labelPlacementMethod="XHeight" showAxis="1" barWidth="5" opacity="1" penAlpha="255" minimumSize="0">
      <fontProperties description="Sans Serif,9,-1,5,50,0,0,0,0,0" style=""/>
      <axisSymbol>
        <symbol clip_to_extent="1" type="line" force_rhr="0" name="" alpha="1">
          <layer enabled="1" class="SimpleLine" pass="0" locked="0">
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
  <DiagramLayerSettings dist="0" priority="0" placement="1" linePlacementFlags="18" obstacle="0" zIndex="0" showAll="1">
    <properties>
      <Option type="Map">
        <Option value="" type="QString" name="name"/>
        <Option name="properties"/>
        <Option value="collection" type="QString" name="type"/>
      </Option>
    </properties>
  </DiagramLayerSettings>
  <geometryOptions geometryPrecision="0" removeDuplicateNodes="0">
    <activeChecks/>
    <checkConfiguration type="Map">
      <Option type="Map" name="QgsGeometryGapCheck">
        <Option value="0" type="double" name="allowedGapsBuffer"/>
        <Option value="false" type="bool" name="allowedGapsEnabled"/>
        <Option value="" type="QString" name="allowedGapsLayer"/>
      </Option>
    </checkConfiguration>
  </geometryOptions>
  <legend type="default-vector"/>
  <referencedLayers/>
  <fieldConfiguration>
    <field name="DOWNWIND_FACADE_ID" configurationFlags="None">
      <editWidget type="Range">
        <config>
          <Option/>
        </config>
      </editWidget>
    </field>
    <field name="ID_STACKED_BLOCK" configurationFlags="None">
      <editWidget type="Range">
        <config>
          <Option/>
        </config>
      </editWidget>
    </field>
    <field name="HEIGHT_ROO" configurationFlags="None">
      <editWidget type="Range">
        <config>
          <Option/>
        </config>
      </editWidget>
    </field>
    <field name="STACKED_BLOCK_X_MED" configurationFlags="None">
      <editWidget type="TextEdit">
        <config>
          <Option/>
        </config>
      </editWidget>
    </field>
    <field name="STACKED_BLOCK_UPSTREAMEST_X" configurationFlags="None">
      <editWidget type="TextEdit">
        <config>
          <Option/>
        </config>
      </editWidget>
    </field>
    <field name="SIN_BLOCK_LEFT_AZIMUTH" configurationFlags="None">
      <editWidget type="TextEdit">
        <config>
          <Option/>
        </config>
      </editWidget>
    </field>
    <field name="COS_BLOCK_LEFT_AZIMUTH" configurationFlags="None">
      <editWidget type="TextEdit">
        <config>
          <Option/>
        </config>
      </editWidget>
    </field>
    <field name="COS_BLOCK_RIGHT_AZIMUTH" configurationFlags="None">
      <editWidget type="TextEdit">
        <config>
          <Option/>
        </config>
      </editWidget>
    </field>
    <field name="SIN_BLOCK_RIGHT_AZIMUTH" configurationFlags="None">
      <editWidget type="TextEdit">
        <config>
          <Option/>
        </config>
      </editWidget>
    </field>
    <field name="BLOCK_WIDTH" configurationFlags="None">
      <editWidget type="TextEdit">
        <config>
          <Option/>
        </config>
      </editWidget>
    </field>
  </fieldConfiguration>
  <aliases>
    <alias name="" field="DOWNWIND_FACADE_ID" index="0"/>
    <alias name="" field="ID_STACKED_BLOCK" index="1"/>
    <alias name="" field="HEIGHT_ROO" index="2"/>
    <alias name="" field="STACKED_BLOCK_X_MED" index="3"/>
    <alias name="" field="STACKED_BLOCK_UPSTREAMEST_X" index="4"/>
    <alias name="" field="SIN_BLOCK_LEFT_AZIMUTH" index="5"/>
    <alias name="" field="COS_BLOCK_LEFT_AZIMUTH" index="6"/>
    <alias name="" field="COS_BLOCK_RIGHT_AZIMUTH" index="7"/>
    <alias name="" field="SIN_BLOCK_RIGHT_AZIMUTH" index="8"/>
    <alias name="" field="BLOCK_WIDTH" index="9"/>
  </aliases>
  <defaults>
    <default field="DOWNWIND_FACADE_ID" applyOnUpdate="0" expression=""/>
    <default field="ID_STACKED_BLOCK" applyOnUpdate="0" expression=""/>
    <default field="HEIGHT_ROO" applyOnUpdate="0" expression=""/>
    <default field="STACKED_BLOCK_X_MED" applyOnUpdate="0" expression=""/>
    <default field="STACKED_BLOCK_UPSTREAMEST_X" applyOnUpdate="0" expression=""/>
    <default field="SIN_BLOCK_LEFT_AZIMUTH" applyOnUpdate="0" expression=""/>
    <default field="COS_BLOCK_LEFT_AZIMUTH" applyOnUpdate="0" expression=""/>
    <default field="COS_BLOCK_RIGHT_AZIMUTH" applyOnUpdate="0" expression=""/>
    <default field="SIN_BLOCK_RIGHT_AZIMUTH" applyOnUpdate="0" expression=""/>
    <default field="BLOCK_WIDTH" applyOnUpdate="0" expression=""/>
  </defaults>
  <constraints>
    <constraint unique_strength="0" notnull_strength="0" field="DOWNWIND_FACADE_ID" exp_strength="0" constraints="0"/>
    <constraint unique_strength="0" notnull_strength="0" field="ID_STACKED_BLOCK" exp_strength="0" constraints="0"/>
    <constraint unique_strength="0" notnull_strength="0" field="HEIGHT_ROO" exp_strength="0" constraints="0"/>
    <constraint unique_strength="0" notnull_strength="0" field="STACKED_BLOCK_X_MED" exp_strength="0" constraints="0"/>
    <constraint unique_strength="0" notnull_strength="0" field="STACKED_BLOCK_UPSTREAMEST_X" exp_strength="0" constraints="0"/>
    <constraint unique_strength="0" notnull_strength="0" field="SIN_BLOCK_LEFT_AZIMUTH" exp_strength="0" constraints="0"/>
    <constraint unique_strength="0" notnull_strength="0" field="COS_BLOCK_LEFT_AZIMUTH" exp_strength="0" constraints="0"/>
    <constraint unique_strength="0" notnull_strength="0" field="COS_BLOCK_RIGHT_AZIMUTH" exp_strength="0" constraints="0"/>
    <constraint unique_strength="0" notnull_strength="0" field="SIN_BLOCK_RIGHT_AZIMUTH" exp_strength="0" constraints="0"/>
    <constraint unique_strength="0" notnull_strength="0" field="BLOCK_WIDTH" exp_strength="0" constraints="0"/>
  </constraints>
  <constraintExpressions>
    <constraint exp="" desc="" field="DOWNWIND_FACADE_ID"/>
    <constraint exp="" desc="" field="ID_STACKED_BLOCK"/>
    <constraint exp="" desc="" field="HEIGHT_ROO"/>
    <constraint exp="" desc="" field="STACKED_BLOCK_X_MED"/>
    <constraint exp="" desc="" field="STACKED_BLOCK_UPSTREAMEST_X"/>
    <constraint exp="" desc="" field="SIN_BLOCK_LEFT_AZIMUTH"/>
    <constraint exp="" desc="" field="COS_BLOCK_LEFT_AZIMUTH"/>
    <constraint exp="" desc="" field="COS_BLOCK_RIGHT_AZIMUTH"/>
    <constraint exp="" desc="" field="SIN_BLOCK_RIGHT_AZIMUTH"/>
    <constraint exp="" desc="" field="BLOCK_WIDTH"/>
  </constraintExpressions>
  <expressionfields/>
  <attributeactions>
    <defaultAction value="{00000000-0000-0000-0000-000000000000}" key="Canvas"/>
  </attributeactions>
  <attributetableconfig sortExpression="" sortOrder="0" actionWidgetStyle="dropDown">
    <columns>
      <column type="field" width="-1" name="DOWNWIND_FACADE_ID" hidden="0"/>
      <column type="field" width="-1" name="ID_STACKED_BLOCK" hidden="0"/>
      <column type="field" width="-1" name="HEIGHT_ROO" hidden="0"/>
      <column type="field" width="-1" name="STACKED_BLOCK_X_MED" hidden="0"/>
      <column type="field" width="-1" name="STACKED_BLOCK_UPSTREAMEST_X" hidden="0"/>
      <column type="field" width="-1" name="SIN_BLOCK_LEFT_AZIMUTH" hidden="0"/>
      <column type="field" width="-1" name="COS_BLOCK_LEFT_AZIMUTH" hidden="0"/>
      <column type="field" width="-1" name="COS_BLOCK_RIGHT_AZIMUTH" hidden="0"/>
      <column type="field" width="-1" name="SIN_BLOCK_RIGHT_AZIMUTH" hidden="0"/>
      <column type="field" width="-1" name="BLOCK_WIDTH" hidden="0"/>
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
    <field editable="1" name="BLOCK_WIDTH"/>
    <field editable="1" name="COS_BLOCK_LEFT_AZIMUTH"/>
    <field editable="1" name="COS_BLOCK_RIGHT_AZIMUTH"/>
    <field editable="1" name="DOWNWIND_FACADE_ID"/>
    <field editable="1" name="HEIGHT_ROO"/>
    <field editable="1" name="ID_STACKED_BLOCK"/>
    <field editable="1" name="SIN_BLOCK_LEFT_AZIMUTH"/>
    <field editable="1" name="SIN_BLOCK_RIGHT_AZIMUTH"/>
    <field editable="1" name="STACKED_BLOCK_UPSTREAMEST_X"/>
    <field editable="1" name="STACKED_BLOCK_X_MED"/>
  </editable>
  <labelOnTop>
    <field labelOnTop="0" name="BLOCK_WIDTH"/>
    <field labelOnTop="0" name="COS_BLOCK_LEFT_AZIMUTH"/>
    <field labelOnTop="0" name="COS_BLOCK_RIGHT_AZIMUTH"/>
    <field labelOnTop="0" name="DOWNWIND_FACADE_ID"/>
    <field labelOnTop="0" name="HEIGHT_ROO"/>
    <field labelOnTop="0" name="ID_STACKED_BLOCK"/>
    <field labelOnTop="0" name="SIN_BLOCK_LEFT_AZIMUTH"/>
    <field labelOnTop="0" name="SIN_BLOCK_RIGHT_AZIMUTH"/>
    <field labelOnTop="0" name="STACKED_BLOCK_UPSTREAMEST_X"/>
    <field labelOnTop="0" name="STACKED_BLOCK_X_MED"/>
  </labelOnTop>
  <dataDefinedFieldProperties/>
  <widgets/>
  <previewExpression>"DOWNWIND_FACADE_ID"</previewExpression>
  <mapTip></mapTip>
  <layerGeometryType>2</layerGeometryType>
</qgis>
