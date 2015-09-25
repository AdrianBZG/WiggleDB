//////////////////////////////////////////
// Global configuration
//////////////////////////////////////////

var CGI_URL = "http://" + location.hostname + "/cgi-bin/wiggleCGI.py?";
var attribute_values_file = "datasets.attribs.json";

//////////////////////////////////////////
// Main function 
//////////////////////////////////////////

$(document).ready(main)

function main() {
  $.getJSON(CGI_URL + "annotations=1").done(get_annotations).fail(catch_JSON_error);
}

function get_annotations(data) {
  annotations = data['annotations'];
  create_all_selectors();
  add_annotations();
  define_buttons();
}

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

var annotations = null;

var user_annotations = null;

var reduction_opts = {"Intersection":"unit mult", "Union":"unit sum"};

var comparison_opts = {"Intersection":"unit mult", "Union":"unit sum", "Difference": "unit diff"};

var annotation_opts = {"Intersection":"unit mult", "Union":"unit sum", "Difference": "unit diff", "Overlap frequency": "overlaps"};

//////////////////////////////////////////
// Creating multiselects 
//////////////////////////////////////////

function add_value_to_multiselect(value, div) {
  $("<option>").attr("value",value).text(value).appendTo(div);
}

function create_multiselect(container, attribute, panel) {
  var multiselect2 = $("<select>")
    .addClass("multiselect")
    .attr("multiple","multiple")
    .appendTo(container)
    .attr("attribute",panel_letters[panel.attr("id")]+ "_" + attribute);

  if (attribute in attribute_values) {
    attribute_values[attribute].sort().map(function(value) {add_value_to_multiselect(value, multiselect2);});
  }
  multiselect2.multiselect({onChange: function(element, checked) {update_panel_count(panel);}, maxHeight: 400, buttonWidth:'100%'});
  multiselect2.parent().find('.btn').css("white-space","normal");
}

function all_selects_are_used(panel) {
  var selects = panel.find(".form-control");
  var res = true; 
  selects.each(function(rank, select) {if ($(select).val() == "None") {res = false;}});
  return res;
}

function change_multiselect() {
  var select = $(this);
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
  Object.keys(attribute_values).sort().map(function(attribute) {add_attribute_to_select(attribute, select);});
  $("<option>").attr("value","None").text("None").attr("selected","selected").appendTo(select);
  select.change(change_multiselect);
}

function add_attribute_to_select(attribute, select) {
  $("<option>").attr("value",attribute).text(attribute).appendTo(select);
}

function remove_attribute_to_select(attribute, select) {
  var search_string = "[value='" + attribute + "']";
  $(select).find(search_string).remove();
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
  var col1 = $("<div>").addClass("form-group").addClass("col-md-5").appendTo(row);
  $("<div>").addClass("form-group").addClass("col-md-1").text(" is: ").appendTo(row);
  var col2 = $("<div>").addClass("form-group").addClass("col-md-5").appendTo(row);

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

  $("<input>").attr('type','text').attr('id','distance').attr('style','width: 50%;').attr('value','0').change(update_my_tab).prependTo(col2);

  var object = $("<select>").addClass("form-control").addClass('filter_reference').attr('id','reference').appendTo(col3);
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
// Creating the personal reference selectors 
//////////////////////////////////////////

function add_personal_annotation_to_references(annotation) {
  var row = $("<tr>").appendTo($('#refs')).attr("value", annotation);
  $("<td>").appendTo(row).append($("<input>").attr("class", "select_annot").attr("type","checkbox").attr("value", annotation));
  $("<td>").appendTo(row).text("\t" + annotation);
  if (user_annotations[annotation]) {
    $("<td>").appendTo(row).attr("align","right").append(
      $("<input>").attr("id","history").attr("type","image").attr("src","images/info.png").attr("value",annotation).attr("title", "Info")
    );
  }
  $("<td>").appendTo(row);
}

function add_personal_annotation_to_personal_list(annotation) {
  var row = $("<tr>").appendTo($('#myannots')).attr("value", annotation);
  $("<td>").appendTo(row).attr("align","right").append(
    $("<input>").attr("id","share").attr("type","image").attr("src","images/share.cropped.png").attr("value",annotation).attr("title", "View, download or share")
  );
  $("<td>").appendTo(row).attr("align","right").append(
    $("<input>").attr("id","trash").attr("type","image").attr("src","images/trash.cropped.png").attr("value",annotation).attr("title", "Delete")
  );
  if (user_annotations[annotation]) {
    $("<td>").appendTo(row).attr("align","right").append(
      $("<input>").attr("id","history").attr("type","image").attr("src","images/info.png").attr("value",annotation).attr("title", "Info")
    );
  } else {
    $("<td>").appendTo(row);
  }

  $("<td>").appendTo(row).text("\t\t" + annotation);
}

function add_personal_annotation(annotation) {
  add_personal_annotation_to_references(annotation)
  add_personal_annotation_to_personal_list(annotation)
}

function remove_personal_annotation(annotation) {
  var search_string = "[value='" + annotation + "']";
  $("#refs").find(search_string).remove();
  $("#myannots").find(search_string).remove();
}

function add_personal_filters(annotation) {
  $(".filter_reference").each(
    function (index, element) {
      add_attribute_to_select(annotation, element);
    }
  );
}

function remove_personal_filters(annotation) {
  $(".filter_reference").each(
    function (index, element) {
      remove_attribute_to_select(annotation, element);
    }
  );
}

function update_default_annotation_name(name, annotation) {
  if (name == annotation) {
    var tail = annotation.substr(14);
    if (tail == "") {
      return "My Annotation 2";
    } else {
      var new_index = parseInt(tail) + 1;
      return "My Annotation " + new_index.toString();
    } 
  } else {
    return name;
  }
}

function update_default_annotation_names(annotation) {
  $("#upload").find("input[id$='description']").each( 
    function (index, element) {
      $(element).attr("value", update_default_annotation_name(name, annotation));
    }
  );
}

function insert_personal_annotation(annotation) {
  add_personal_annotation(annotation);
  add_personal_filters(annotation);
  update_default_annotation_names(annotation);
}

function insert_personal_annotations(data) {
  user_annotations = data['files'];
  Object.keys(user_annotations).map(insert_personal_annotation);
}

function add_personal_annotations() {
  $("#myannots").on("click","#trash", confirm_remove_personal_file);
  $("#myannots").on("click","#share", share_personal_file);
  $("#myannots").on("click","#history", history);
  $("#refs").on("click","#history", history);
  jQuery.getJSON(CGI_URL + "myannotations=1&userid=" + user_ID()).done(insert_personal_annotations).fail(catch_JSON_error);
}

function personal_annotation_deleted(annotation) {
  remove_personal_annotation(annotation);
  remove_personal_filters(annotation);
}

function remove_personal_file(annotation) {
  jQuery.getJSON(CGI_URL + "&userid=" + user_ID() + "&remove_annotation=" + annotation).done(function (data) {personal_annotation_deleted(annotation)}).fail(catch_JSON_error);
}

function confirm_remove_personal_file() {
  var annotation = $(this).attr('value');
  var modal = $('#Confirm_deletion_modal').clone();
  modal.find('#name').text(annotation);
  modal.find('#confirm').on("click", function () {remove_personal_file(annotation)})
  modal.modal();
}

function display_annotation_location(data) {
  var modal = $('#Share_modal').clone();
  modal.find('#name').text(data['name']);
  modal.find('#url').text(data['url']);
  modal.find('#url').attr('href', data['url']);
  modal.modal();
}

function share_personal_file() {
  var annotation = $(this).attr('value');
  jQuery.getJSON(CGI_URL + "&userid=" + user_ID() + "&share_annotation=" + annotation).done(display_annotation_location).fail(catch_JSON_error);
}

function display_history(data) {
  var modal = $('#History_modal').clone();
  modal.find('#myModalLabel').text(data['name']);
  modal.find('#name').text(data['name']);
  modal.find('#history').text(data['history']);
  modal.modal();
}

function history() {
  var annotation = $(this).attr('value');
  jQuery.getJSON(CGI_URL + "&userid=" + user_ID() + "&history=" + annotation).done(display_history).fail(catch_JSON_error);
}

//////////////////////////////////////////
// Creating the reference selectors 
//////////////////////////////////////////

function display_provenance(data) {
  var modal = $('#Provenance_modal').clone();
  modal.find('#myModalLabel').text(data["name"]);
  modal.find('#myModalBody').text(data["description"]);
  modal.modal();
}

function get_provenance() {
  var dataset = $(this).attr('value');
  jQuery.getJSON(CGI_URL + "provenance=" + dataset).done(display_provenance).fail(catch_JSON_error);
}

function add_annotation(annotation) {
  var row = $("<tr>").appendTo($('#refs'));
  $("<td>").appendTo(row).append($("<input>").attr("class", "select_annot").attr("type","checkbox").attr("value", annotation));
  $("<td>").appendTo(row).text("\t" + annotation);
  $("<td>").appendTo(row).attr("align","right").append($("<input>").attr("id","source").attr("type","image").attr("src","images/info.png").attr("value",annotation));
}

function add_annotations() {
  $("#refs").on("click","#source", get_provenance);
  fill_select($('#reference_reduction'), reduction_opts);
  annotations.map(add_annotation);
  if (getCookie("username") != "" && getCookie("username") != null) {
    add_personal_annotations();
  }
}

///////////////////////////////////////////
// Panel filters
///////////////////////////////////////////

function panel_filters(panel) {
  var filters = [];
  var header = "filter_" + panel_letters[panel.attr("id")] + "=";
  panel.find('.filter').each(
    function (index, div) {
      var constraint = $(div).find("#constraint").val();
      var distance = $(div).find("#distance").val();
      var reference = $(div).find("#reference").val();
      if (constraint == "(No filter)") {
	return;
      } else if (distance == null || distance == "") {
	filters.push(header + [constraint, "0", reference].join("|")); 
      } else {
	filters.push(header + [constraint, distance, reference].join("|")); 
      } 
    }
  );
  return filters.join("&");
}

//////////////////////////////////////////
// Panel Query
//////////////////////////////////////////

function panel_query(panel) {
  var string = [];
  panel.find(".multiselect option:selected").each(
    function (index, element) {
      var variable = $(element).parent().attr("attribute");
      var value = $(element).val();  
      if (value != null && value != "") {
        string.push(variable+'='+value);
      } 
    }
  );
  return string.join("&");
}

///////////////////////////////////////////
// Panel reduction
///////////////////////////////////////////

function panel_reduction(panel) {
  return panel.find('#reduction').val();
}

//////////////////////////////////////////
// Computing panel selection count
//////////////////////////////////////////

function update_panel_count(panel) {
  $.getJSON(CGI_URL + "count=true&" + panel_query(panel)).done(
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

function update_panel_reduction(panel) {
  fill_select(panel.find('#reduction'), reduction_opts)
}

//////////////////////////////////////////
// Updating tab comparison select 
//////////////////////////////////////////

function update_tab_comparison(tab) {
  if (tab.attr('id') == 'annotate') {
    fill_select(tab.find("#annotation"), annotation_opts);
  } else {
    fill_select(tab.find("#comparison"), comparison_opts);
  }
}

//////////////////////////////////////////
// Button up!
//////////////////////////////////////////

function define_buttons() {
  $('#summary_button').click(summary);
  $('#comparison_button').click(comparison);
  $('#annotation_button').click(annotation);
  $('#upload_button').click(upload_url);
  $('#upload_file_button').click(upload_file);
  reset_buttons();
}

function reset_buttons() {
  $('#summary_button').button('reset');
  $('#comparison_button').button('reset');
  $('#annotation_button').button('reset');
  $('#upload_button').button('reset');
  $('#upload_file_button').button('reset');
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

function save_result_file() {
  var modal = $(this).parents('.modal');
  var name = modal.find('#name').val();
  var url = modal.find('#url').attr("href");
  user_annotations[name] = 1;
  $.getJSON(CGI_URL + "userid=" + user_ID() + "&uploadUrl=" + url + "&description=" + name).done(report_upload).fail(catch_JSON_error);
}

function report_result(data) {
  if (data["status"] == "DONE") {
    if (data['url'].substr(-4,4) == ".txt") {
      var modal = $('#Image_modal').clone();
      modal.find('#url').attr('href',data['url']);
      modal.find('img').attr('src',data['view']);
      modal.find('#photo_url').attr('href',data['view']);
      modal.modal();
    } else if (data['url'].substr(-3,3) == '.bw' || data['url'].substr(-3,3) == '.bb' || data['url'].substr(-4,4) == '.bed') {
      var modal = $('#Success_modal').clone();
      modal.find('#url').attr('href',data['url']);
      modal.find('#view').attr('href',data['view']);
      modal.on('click','#save', save_result_file);
      modal.modal();
    }
  } else if (data["status"] == "EMPTY") {
    $('#Empty_modal').modal();	
  } else if (data["status"] == "INVALID") {
    $('#Invalid_modal').modal();	
  } else if (data["status"] == "ERROR") {
    $('#Failure_modal').modal();	
  }
  reset_buttons();
}

// Get result
function get_result() {
  $.getJSON(CGI_URL + "result=" + $('#result_box').val()).done(report_result).fail(catch_JSON_error);
}

///////////////////////////////////////////
// Button actions
///////////////////////////////////////////

function upload_file() {
  $(this).button('loading');
  var formData = new FormData();
  formData.append("file", $(this).parents("#FileUpload").find("#file").get(0).files[0]);
  $.ajax({
    url: CGI_URL  + "uploadFile=1&userid=" + user_ID() + "&description=" + $(this).parents("#FileUpload").find("#upload_description").val(),
    data: formData, 
    processData: false, 
    type: 'POST',
    contentType: false,
    success: report_upload,
    fail: catch_JSON_error
  });
  
}

// Upload user dataset
function upload_url() {
  var files = $("#upload").find("#file").val();
  $(this).button('loading');
  $.getJSON(CGI_URL + "userid=" + user_ID() + "&uploadUrl=" + $('#URL').val() + "&description=" + $("#url_description").val()).done(report_upload).fail(catch_JSON_error);
}

function report_upload(data) {
  reset_buttons();
  if (data["status"] == "UPLOADED") {
      user_annotations[data['name']] = data['has_history'];
      insert_personal_annotation(data['name'])
      var modal = $('#Success_upload_modal').clone();
      modal.modal();
  } else if (data['status'] == 'NAME_USED') {
      var modal = $('#Rename_modal').clone();
      modal.find('#url').text(data['url']);
      modal.modal();
  } else if (data['status'] == 'WRONG_ASSEMBLY') {
      var modal = $('#Assembly_modal').clone();
      modal.find('#expected_chromosomes').text(data['expected_chromosomes'].join(", "));
      modal.find('#found_chromosomes').text(data['found_chromosomes'].join(", "));
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
  $.getJSON(CGI_URL + query).done(report_result).fail(catch_JSON_error);
}

function get_optional_user_id() {
  var user = getCookie("username");
  if (user != "") {
    return "userid=" + user;
  } else {
    return "";
  }
}

// Request summary
function summary() {
  var panel = $('#choose');
  $(this).button('loading');
  submit_query([panel_query(panel), panel_filters(panel), 'wa=' + panel_reduction(panel), get_optional_user_id()].join("&")); 
}

// Request comparison
function comparison() {
  $(this).button('loading');
  var comparison = $('#comparison').val();
  var panelA = $('#chooseA');
  var panelB = $('#chooseB');
  submit_query(
    [ 
      panel_query(panelA),
      panel_filters(panelA),
      panel_query(panelB),
      panel_filters(panelB),
      get_optional_user_id(),
      'wa='+panel_reduction(panelA),
      'wb='+panel_reduction(panelB),
      'w='+comparison
    ].join("&")
  ); 
}

// Request annotation
function annotation() {
  $(this).button('loading');
  var comparison = $('#annotation').val();
  var panelA = $('#chooseA2');
  var annots = [];
  $('#refs').find(".select_annot").each(function (rank, obj) {if (obj.checked) annots.push($(obj).attr("value"));});
  submit_query(
    [ 
      panel_query(panelA),
      annots.map(function (str) {return "B_annot_name="+str}).join("&"),
      get_optional_user_id(),
      'wa='+panel_reduction(panelA),
      'w='+comparison
    ].join("&")
  ); 
}

// JSON error handler
function catch_JSON_error(jqXHR, textStatus, errorThrown) {
  console.log('JSON failed: ' + textStatus + ":" + errorThrown + ":" + jqXHR.responseText + ":" + jqXHR.getAllResponseHeaders());
}

///////////////////////////////////////////
// The Cookie Jar
///////////////////////////////////////////

function setCookie(cname, cvalue, exdays) {
  var d = new Date();
  d.setTime(d.getTime() + (exdays*24*60*60*1000));
  var expires = "expires="+d.toUTCString();
  document.cookie = cname + "=" + cvalue + "; " + expires;
}

function getCookie(cname) {
  var name = cname + "=";
  var ca = document.cookie.split(';');
  for(var i=0; i<ca.length; i++) {
    var c = ca[i];
    while (c.charAt(0)==' ') c = c.substring(1);
    if (c.indexOf(name) == 0) return c.substring(name.length,c.length);
  }
  return "";
} 

///////////////////////////////////////////
// User ID
///////////////////////////////////////////

function generate_ID() {
  return 'xxxxxxxx-xxxx-4xxx-yxxx-xxxxxxxxxxxx'.replace(
    /[xy]/g, 
    function(c) {
      var r = Math.random()*16|0, v = c == 'x' ? r : (r&0x3|0x8);
      return v.toString(16);
    }
  );
}

function user_ID() {
  var user = getCookie("username");
  if (user != "") {
    return user;
  } else {
    setCookie("username", generate_ID(), 180);
  }    
}
