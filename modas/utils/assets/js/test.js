$(function() {
  function clickNav(event) {
    var target = event.currentTarget;
    var tab = $(target).attr("href");
    // console.log(tab);
    $('a[href="' + tab + '"]').tab('show');
    $(tab).tab('show');
    if ($("a[role='tab']")) {
      var arrLink = $("a[role='tab']");
      for (i = 0; i < arrLink.length; i++) {
        var link = arrLink[i];
        if (link.classList.contains('active')) {
          link.classList.remove('active');
        }
      }
    }
    target.classList.add('active');
  }

  $(".collapse").on('show.bs.collapse', function() {
    $(this).prev(".qtlItem").find(".fa").removeClass("fa-plus").addClass("fa-minus");
  }).on('hide.bs.collapse', function() {
    $(this).prev(".qtlItem").find(".fa").removeClass("fa-minus").addClass("fa-plus");
    $(this).prev(".qtlItem").removeClass("active");
  });

  if ($("a[role='tab']")) {
    var arrLink = $("a[role='tab']");
    for (i = 0; i < arrLink.length; i++) {
      var link = arrLink[i];
      link.classList.remove('active');
      link.onclick = clickNav;
    }
  }
});
