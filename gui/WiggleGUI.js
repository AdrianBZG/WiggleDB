//////////////////////////////////////////
// Global configuration
//////////////////////////////////////////

var CGI_URL = "http://" + location.hostname + "/cgi-bin/wiggleCGI.py?";
var attribute_values_file = "datasets.attribs.json";
var assembly = "GRCh37";

//////////////////////////////////////////
// Main function 
//////////////////////////////////////////

$(document).ready(
  function () {
    create_all_selectors();
    add_annotations();
    define_buttons();
  }
)

//////////////////////////////////////////
// Global variables
//////////////////////////////////////////

var selection_panels = [
  "choose",
  "chooseB",
  "chooseA2",
  "chooseA"
];

var panel_letters = {
  "choose":"A",
  "chooseB":"B",
  "chooseA2":"A",
  "chooseA":"A"
};

var attribute_values = null;

var reduction_opts = {
  "signal": {"Sum":"sum","Mininum":"min","Maximum":"max","Mean":"mean","Median":"median"},
  "regions": {"Intersection":"unit mult", "Union":"unit sum"}
};

var comparison_opts = {
  "regions": {"Intersection":"unit mult", "Union":"unit sum", "Difference": "unit diff"},
  "signal": {"Difference":"diff", "Ratio":"ratio", "Log Ratio":"ln ratio", "T-test":"t-test", "Wilcoxon rank test":"wilcoxon"},
  "mixed": {"Distribution":"histogram", "Profile curve":"profile", "Profile matrix":"profiles"}
};

var annotation_opts = {
  "regions": {"Intersection":"unit mult", "Union":"unit sum", "Difference": "unit diff", "Overlap frequency": "overlaps"},
  "signal": {"Distribution":"histogram", "Profile curve":"profile", "Profile matrix":"profiles", "Enrichment":"apply_paste"}
};

//////////////////////////////////////////
// Creating multiselects 
//////////////////////////////////////////

function add_value_to_multiselect(value, div) {
  $("<option>").attr("value",value).text(value).appendTo(div);
}

function all_selects_are_used(panel) {
  var selects = panel.find(".form-control");
  var res = true; 
  selects.each(function(rank, select) {if ($(select).val() == "None") {res = false;}});
  return res;
}

function create_multiselect(container, attribute, panel) {
  var multiselect2 = $("<select>")
    .addClass("multiselect")
    .attr("multiple","multiple")
    .appendTo(container)
    .attr("attribute",panel_letters[panel.attr("id")]+ "_" + attribute);

  if (attribute in attribute_values) {
    attribute_values[attribute].map(function(value) {add_value_to_multiselect(value, multiselect2);});
  }
  multiselect2.multiselect({onChange: function(element, checked) {update_panel_count(panel);}, maxHeight: 400, buttonWidth:'100%'});
  multiselect2.parent().find('.btn').css("white-space","normal");
}

function reset_multiselect(select) {
  var panel = select.parents("[id*='choose']");
  var col = select.parent().parent().find(".multiselect").parent(".form-group");
  if (select.val() in attribute_values) {
    col.children().remove();
    create_multiselect(col, select.val(), panel);
  } else {
    col.parent().remove();
  }
  
  if (all_selects_are_used(panel)) {
    create_selection_div(panel);
  }
}

function create_attribute_select(container) {
  var select = $("<select>").addClass("form-control").appendTo(container);
  Object.keys(attribute_values).map(function(attribute) {add_attribute_to_select(attribute, select);});
  $("<option>").attr("value","None").text("None").attr("selected","selected").appendTo(select);
  select.change(function() {reset_multiselect(select);})
}

function add_attribute_to_select(attribute, select) {
  $("<option>").attr("value",attribute).text(attribute).appendTo(select);
}

function update_panel(panel) {
  // Compute initial count:
  update_panel_count(panel);
  // Set reduction select options
  update_panel_reduction(panel);
  // Comparison select:
  update_tab_comparison(panel.parents('.tab-pane'));
}

function create_selection_div(panel) {
  // Top most container
  var row = $("<div>").addClass("row").appendTo(panel.find("#selection"));

  // Division into fixed width columns
  var col1 = $("<div>").addClass("form-group").addClass("col-lg-5").appendTo(row);
  $("<div>").addClass("form-group").addClass("col-lg-1").text(" is: ").appendTo(row);
  var col2 = $("<div>").addClass("form-group").addClass("col-lg-5").appendTo(row);

  // Create attribute selector in column 1:
  create_attribute_select(col1);
  // Create empty value selector in column 2:
  create_multiselect(col2, "None", panel);

  update_panel(panel);
}

function all_filters_used(panel) {
  var res = true;
  panel.find(".filter").each(
    function (index, obj) {
      if ($(obj).find("#constraint").val() == "(No filter)") {
        res = false;
      }
    }
  );
  return res;
}

function reset_filter() {
  var panel = $(this).parents(".tab-pane");
  if ($(this).val() == "(No filter)" && panel.find(".filter").length > 1) {
    $(this).parents(".filter").remove();
  }
  if (all_filters_used(panel)) {
    create_filter_div(panel);
  }
}

function create_filter_div(panel) {
  var row = $('<div>').attr('class','filter row').appendTo(panel.find(".filters"));

  // Division into fixed width columns
  var col1 = $("<div>").addClass("form-group").addClass("col-md-4").appendTo(row);
  var col2 = $("<div>").addClass("form-group").text(" bp from ").addClass("col-md-4").appendTo(row);
  var col3 = $("<div>").addClass("form-group").addClass("col-md-4").appendTo(row);

  var verb = $("<select>").addClass("form-control").attr('id','constraint').appendTo(col1);
  $("<option>").attr("value","overlaps").text("are within").attr("selected","selected").appendTo(verb);
  $("<option>").attr("value","noverlaps").text("are farther than").attr("selected","selected").appendTo(verb);
  $("<option>").attr("value",null).text("(No filter)").attr("selected","selected").appendTo(verb);

  $("<input>").attr('type','text').attr('id','distance').attr('style','width: 50%;').change(update_my_tab).prependTo(col2);

  var object = $("<select>").addClass("form-control").attr('id','reference').appendTo(col3);
  if (annotations != null) {
    annotations.map(function(attribute) {add_attribute_to_select(attribute, object);});
  }

  verb.change(reset_filter);
}

function create_panel(panel) {
  create_selection_div(panel);
  create_filter_div(panel);
}

function create_all_selectors() {
  jQuery.getJSON(attribute_values_file).done(function(values) {
    attribute_values = values;
    selection_panels.map(function (id) {create_panel($("#"+id));});
  }).fail(catch_JSON_error);
}

//////////////////////////////////////////
// Creating the reference selectors 
//////////////////////////////////////////

function provenance(data) {
  var modal = $('#Provenance_modal').clone();
  modal.find('#myModalLabel').text(data["name"]);
  modal.find('#myModalBody').text(data["description"]);
  modal.modal();
}

function get_provenance() {
  var dataset = $(this).attr('value');
  jQuery.getJSON(CGI_URL + "provenance=" + dataset).done(provenance).fail(catch_JSON_error);
}

function add_annotation(annotation) {
  var row = $("<tr>").appendTo($('#refs'));
  $("<td>").appendTo(row).append($("<input>").attr("class", "select_annot").attr("type","checkbox").attr("value", annotation));
  $("<td>").appendTo(row).text("\t" + annotation);
  $("<td>").appendTo(row).attr("align","right").append($("<input>").attr("id","source").attr("type","image").attr("src","images/info.png").attr("value",annotation));
}

function add_annotations_2(data) {
  data['annotations'].map(add_annotation);
}

function add_annotations() {
  //$.getJSON(CGI_URL + "assembly=" + assembly + "&annotations=1").done(add_annotations_2).fail(catch_JSON_error);
  $("#refs").on("click","#source", get_provenance);
  fill_select($('#reference_reduction'), reduction_opts['regions']);
  $.getJSON(CGI_URL + "annotations=1").done(add_annotations_2).fail(catch_JSON_error);
}

//////////////////////////////////////////
// Panel Query
//////////////////////////////////////////

function add_multi_select_query(list, multiselect) {
  var variable = multiselect.attr("attribute");
  var values = multiselect.val();  
  if (values != null && values != "") {
    values.map(function(value) {list.push([variable+'='+value])});
  }
  return list;
}

function  get_panel_type(panel) {
  return "regions";
  //return panel.find('#type').find('.active').find('input').attr('value');
}

function get_panel_letter(panel) {
  return panel_letters[panel.attr('id')];
}

function panel_query(panel) {
  list = []
  panel.find(".multiselect").each(function(i,value) {add_multi_select_query(list, $(value));});
  return list.join("&") + "&" + get_panel_letter(panel) + "_type=" + get_panel_type(panel);
}

///////////////////////////////////////////
// Panel reduction
///////////////////////////////////////////

function panel_reduction(panel) {
  var commands = [];
  panel.find('.filter').each(
    function (index, div) {
      var constraint = $(div).find("#constraint").val();
      var distance = $(div).find("#distance").val();
      var reference = $(div).find("#reference").val();
      if (constraint == "(No filter)") {
	return;
      } else if (distance == null || distance == "") {
	commands.push([constraint, reference].join(" ")); 
      } else {
	commands.push([constraint, "extend", distance, reference].join(" ")); 
      } 
    }
  );
  commands.push(panel.find('#reduction').val());
  return commands.join(" ");
}

//////////////////////////////////////////
// Computing panel selection count
//////////////////////////////////////////

function update_panel_count(panel) {
  //url = CGI_URL + "count=true&assembly=" + assembly + "&" + panel_query(panel);
  url = CGI_URL + "count=true&" + panel_query(panel);
  $.getJSON(url).done(
   function(data, textStatus, jqXHR) {
     panel.find("#count").text("(" + data["count"] + " elements selected)");
   }
  ).fail(catch_JSON_error);
}

//////////////////////////////////////////
// Updating panel reduction select 
//////////////////////////////////////////

function fill_select(select, options) {
  select.children().remove();
  Object.keys(options).map(function(opt) {$("<option>").attr("value",options[opt]).text(opt).appendTo(select);});
}

function update_my_tab() {
  update_tab_comparison($(this).parents('.tab-pane'));
}

function add_threshold_selector(div) {
  if (div.find('#threshold').length == 0) {
    var threshold = $('<div>').attr('id','threshold').text("Optional: Select regions above threshold : ").appendTo(div);
    $('<br>').prependTo(threshold);
    $("<input>").attr('type','text').attr('id','threshold_val').change(update_my_tab).appendTo(threshold);
  }
}

function update_panel_reduction(panel) {
  var type = get_panel_type(panel);
  var reduction = panel.find('#reduction');
  fill_select(reduction, reduction_opts[type])
  if (type == 'signal') {
    add_threshold_selector(reduction.parent());
  } else {
    reduction.parent().find('#threshold').remove();
  }
}

//////////////////////////////////////////
// Updating tab comparison select 
//////////////////////////////////////////

function add_width_selector(div) {
  var width = $('<div>').attr('id','width').text("Resolution (points) : ").appendTo(div);
  $('<br>').prependTo(width);
  $("<input>").attr('type','text').attr('id','width_val').attr('value','100').appendTo(width);
}

function update_tab_annotation(tab) {
  var count = 0;  
  tab.find("[id*=choose]").each(function (index, panel) {
    if (get_panel_type($(panel)) == "regions" 
        || $(panel).find('#threshold_val').val() != "") {
      count += 1;
    }
  });
  var select = tab.find("#annotation");
  var div = select.parent().parent();

  if (count == 1) {
    fill_select(select, annotation_opts["regions"]);
    div.find("#width").remove();
  } else {
    fill_select(select, annotation_opts["signal"]);
    div.find("#width").remove();
    add_threshold_selector(div);
  }
}

function update_tab_comparison(tab) {
  if (tab.attr('id') == 'annotate') {
    update_tab_annotation(tab);
  }
  if (tab.attr('id') != 'compare') {
    return;
  } 
  var count = 0;  
  tab.find("[id*=choose]").each(function (index, panel) {
    if (get_panel_type($(panel)) == "regions" 
        || $(panel).find('#threshold_val').val() != "") {
      count += 1;
    }
  });
  var select = tab.find("#comparison");

  // Sneaky way of detecting that we are in annotation tab, not comparison tab
  if (select.length == 0) {
    select = tab.find("#annotation");
    count += 1;
  }

  var div = select.parent().parent();

  if (count == 2) {
    fill_select(select, comparison_opts["regions"]);
    div.find("#width").remove();
    div.find("#threshold").remove();
  } else if (count == 1) {
    fill_select(select, comparison_opts["mixed"]);
    if (div.find("#width").length == 0) {
      add_width_selector(div);
    }
    div.find("#threshold").remove();
  } else {
    fill_select(select, comparison_opts["signal"]);
    div.find("#width").remove();
    add_threshold_selector(div);
  }
  
}

//////////////////////////////////////////
// Update panels
//////////////////////////////////////////

function update_panels() {
  var val = $(this).val();
  if (val != 't-test' && val != 'wilcoxon') {
    $(this).parents('.tab-pane').find('#reduction').prop('disabled', false);
  } else {
    $(this).parents('.tab-pane').find('#reduction').prop('disabled', true);
  }
}

//////////////////////////////////////////
// Button up!
//////////////////////////////////////////

function define_buttons() {
  $('#annotation').change(update_panels);
  $('#comparison').change(update_panels);
  $('#summary_button').click(summary);
  $('#comparison_button').click(comparison);
  $('#annotation_button').click(annotation);
  //$('#result_button').click(get_result);
  $('#upload_button').click(upload_dataset);
  $('[id*=type]').find('label').change(update_my_panel);
}

function update_my_panel() {
  $(this).addClass('active');
  var panel = $(this).parents('[id*="choose"]');
  if (panel.length == 0) {
    return;
  }
  var X = 1;
  update_panel(panel);
  $(this).removeClass('active');
}

function report_result(data) {
  if (data["status"] == "DONE") {
    if (data['url'].substr(-4,4) == ".txt") {
      var modal = $('#Image_modal').clone();
      modal.find('#url').attr('href',data['url']);
      modal.find('img').attr('src',data['view']);
      modal.find('#photo_url').attr('href',data['view']);
      modal.modal();
    } else if (data['url'].substr(-3,3) == '.bw' || data['url'].substr(-4,4) == '.bb') {
      var modal = $('#Success_modal').clone();
      modal.find('#url').attr('href',data['url']);
      modal.find('#view').attr('href',data['view']);
      modal.modal();
    }
  } else if (data["status"] == "EMPTY") {
    $('#Empty_modal').modal();	
  } else if (data["status"] == "INVALID") {
    $('#Invalid_modal').modal();	
  } else if (data["status"] == "ERROR") {
    $('#Failure_modal').modal();	
  } else if (data["status"] == "UNKNOWN") {
    $('#Unknown_modal').modal();	
  } else if (data['status'] == "LAUNCHED") {
    var modal = $("#JobSent_modal").clone();
    modal.find("#job_id").text(data["ID"]);
    modal.modal();
  } else {
    $('#Waiting_modal').modal();	
  }
}

// Get result
function get_result() {
  $.getJSON(CGI_URL + "result=" + $('#result_box').val()).done(report_result).fail(catch_JSON_error);
}

// Upload user dataset
function upload_dataset() {
  $.getJSON(CGI_URL + "uploadUrl=" + $('#uploadURL').val() + "&description=" + $("#uploadDescription").val()).done(report_upload).fail(catch_JSON_error);
}

function report_upload(data) {
  if (data["status"] == "UPLOADED") {
      var modal = $('#Success_upload_modal').clone();
      modal.modal();
  } else if ("format" in data) {
      var modal = $('#Malformed_upload_modal').clone();
      modal.find('#url').text(data['url']);
      modal.find('#format').text(data['format']);
      modal.modal();
  } else {
      var modal = $('#Failed_upload_modal').clone();
      modal.find('#url').text(data['url']);
      modal.find('#format').text(data['format']);
      modal.modal();
  }
}

// Send job to server 
function submit_query(query) {
  //$.getJSON(CGI_URL + "assembly=" + assembly + "&" + query).done(report_result).fail(catch_JSON_error);
  $.getJSON(CGI_URL + query).done(report_result).fail(catch_JSON_error);
}

// Get list of emails separated by ampersands
function emails(panel) {
  return "";
  //var array = panel.find('#email').val().split(/[ ;,]+/);
  //return array.map(function (value, index, array) {return "email="+value;}).join("&");
}

// Request summary
function summary() {
  var panel = $('#choose');
  submit_query(panel_query(panel) + '&wa=' + panel_reduction(panel) + "&" + emails(panel)); 
}

// Request comparison
function comparison() {
  var comparison = $('#comparison').val();
  var width = $('#comparison').parent().parent().find('#width_val');
  if (width.length) {
    comparison = comparison + " " + width.val();
  }
  var threshold = $('#comparison').parent().parent().find('#threshold_val');
  if (threshold.length && threshold.val() != "") {
    comparison = "gt " + threshold.val() + " " + comparison;
  }
  var panelA = $('#chooseA');
  var panelB = $('#chooseB');
  submit_query(
    [ 
      panel_query(panelA),
      panel_query(panelB),
      'wa='+panel_reduction(panelA),
      'wb='+panel_reduction(panelB),
      'w='+comparison,
      emails($('#compare'))
    ].join("&")
  ); 
}

// Request annotation
function annotation() {
  var comparison = $('#annotation').val();
  var width = $('#annotation').parent().parent().find('#width_val')
  if (width.length) {
    comparison = comparison + " " + width.val()
  }
  var panelA = $('#chooseA2');
  var annots = [];
  $('#refs').find(".select_annot").each(function (rank, obj) {if (obj.checked) annots.push($(obj).attr("value"));});
  submit_query(
    [ 
      panel_query(panelA),
      annots.map(function (str) {return "B_annot_name="+str}).join("&"),
      'B_type=regions',
      'wa='+panel_reduction(panelA),
      'w='+comparison,
      emails($('#annotate'))
    ].join("&")
  ); 
}

// JSON error handler
function catch_JSON_error(jqXHR, textStatus, errorThrown) {
  console.log('JSON failed: ' + textStatus + ":" + errorThrown + ":" + jqXHR.responseText + ":" + jqXHR.getAllResponseHeaders());
}
